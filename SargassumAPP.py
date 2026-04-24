import streamlit as st
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import calendar
import requests
from PIL import Image
import io
import cbsodata

# -----------------------------------------------------------------------------
# 1. INITIALIZATION & SETUP
# -----------------------------------------------------------------------------
st.set_page_config(page_title="Sargassum Server-Side Processing App", layout="wide")
st.title("🌊 Sargassum Early Warning & Tourism Regression")
st.markdown("""
**Zero-Download Architecture**: Uses the Copernicus Data Space Ecosystem (CDSE) API.
- **Smart Date Picker**: Queries the CDSE Catalogue to only show dates where satellite data exists.
- **Thresholds Restored**: Dynamically adjust FAI/MCI thresholds inside the server-side evalscript.
- **Robust CBS Parsing**: Dynamically adapts to Statistics Netherlands API changes.
""")

# -----------------------------------------------------------------------------
# 2. DATA COLLECTION: CBS OPEN DATA API
# -----------------------------------------------------------------------------
def safe_parse_cbs(table_id, filter_bonaire=True):
    """
    Safely parses CBS tables regardless of column name changes.
    CBS uses Dutch column names (e.g. 'Eilanden', 'Perioden', 'AantalAankomsten_1')
    even for tables with 'ENG' suffixes, so we cast a wide net.
    """
    try:
        df = pd.DataFrame(cbsodata.get_data(table_id))
        if df.empty:
            return pd.DataFrame()

        # Debug: expose actual columns so users can diagnose mismatches
        with st.expander(f"🔍 Debug: CBS table `{table_id}` raw columns", expanded=False):
            st.write("Columns:", list(df.columns))
            st.write("Sample row:", df.iloc[0].to_dict() if len(df) > 0 else "empty")

        # 1. Find Region/Island column — CBS uses Dutch 'Eilanden' or English 'Islands'
        if filter_bonaire:
            region_col = next(
                (c for c in df.columns if any(x in c.lower() for x in
                 ['territory', 'island', 'eiland', 'region', 'caribbean', 'gebied'])),
                None
            )
            if region_col:
                df = df[df[region_col].astype(str).str.contains('Bonaire', case=False, na=False)]
            # If no region column found, keep all rows (single-island table)
            if df.empty:
                st.warning(f"Table {table_id}: region filter removed all rows. "
                           f"Column tried: {region_col}. Check debug expander above.")
                return pd.DataFrame()

        # 2. Find Period/Date column — CBS uses 'Perioden' (Dutch) or 'Periods' (English)
        period_col = next(
            (c for c in df.columns if any(x in c.lower() for x in ['period', 'perioden', 'datum', 'date'])),
            None
        )
        if not period_col:
            st.warning(f"Table {table_id}: could not find a period/date column. "
                       f"Available columns: {list(df.columns)}")
            return pd.DataFrame()

        # Filter to monthly rows only (CBS monthly format: "2023MM04")
        df = df[df[period_col].astype(str).str.contains('MM', na=False)].copy()
        df['Date'] = pd.to_datetime(
            df[period_col].astype(str).str.replace('MM', '', regex=False),
            format='%Y%m',
            errors='coerce'
        )
        df = df.dropna(subset=['Date'])
        if df.empty:
            st.warning(f"Table {table_id}: no monthly rows found in column '{period_col}'. "
                       "Sample values: " + str(df[period_col].head(5).tolist() if period_col in df.columns else "N/A"))
            return pd.DataFrame()

        # 3. Find value/visitor column
        # CBS Dutch: AantalAankomsten (arrivals), Passagiers (passengers), Aantal (count)
        # CBS English variants: Arrivals, Visitors, Passengers, Total
        candidate_cols = [c for c in df.columns if c not in ['Date', 'ID', period_col]]
        val_col = next(
            (c for c in candidate_cols if any(x in c.lower() for x in
             ['visitor', 'passenger', 'passagier', 'value', 'arrival', 'aankomst',
              'aantal', 'arriv', 'total', 'totaal'])),
            None
        )
        if val_col:
            df['Value'] = pd.to_numeric(df[val_col], errors='coerce').fillna(0)
        else:
            # Last resort: first numeric column that isn't ID
            num_cols = df[candidate_cols].select_dtypes(include=np.number).columns.tolist()
            if num_cols:
                st.info(f"Table {table_id}: guessing value column as '{num_cols[0]}' "
                        f"(from candidates: {candidate_cols}). Check debug expander to verify.")
                df['Value'] = pd.to_numeric(df[num_cols[0]], errors='coerce').fillna(0)
            else:
                st.warning(f"Table {table_id}: no numeric value column found. "
                           f"Columns available: {list(df.columns)}")
                return pd.DataFrame()

        return df[['Date', 'Value']].set_index('Date')

    except Exception as e:
        st.warning(f"Failed to parse CBS table {table_id}: {e}")
        return pd.DataFrame()


@st.cache_data(show_spinner="Loading tourism data from Statistics Netherlands...")
def load_tourism_data():
    """Fetches inbound Air and Cruise tourism safely."""
    df_air = safe_parse_cbs('83104ENG')

    # BUG FIX: guard against empty df_air before renaming
    if df_air.empty:
        st.warning("Air visitor data (table 83104ENG) could not be loaded.")
        return pd.DataFrame()

    df_air = df_air.rename(columns={'Value': 'Air_Visitors'})

    # Find cruise table dynamically
    try:
        tables = pd.DataFrame(cbsodata.get_table_list())
        cruise_tables = tables[tables['Title'].str.contains('Cruise passengers', case=False, na=False)]
    except Exception as e:
        st.warning(f"Could not retrieve CBS table list: {e}")
        cruise_tables = pd.DataFrame()

    if not cruise_tables.empty:
        df_cruise = safe_parse_cbs(cruise_tables.iloc[0]['Identifier'])
        if not df_cruise.empty:
            df_cruise = df_cruise.rename(columns={'Value': 'Cruise_Visitors'})
            df_combined = df_air.join(df_cruise, how='outer').fillna(0)
            df_combined['Total_Bonaire_Visitors'] = df_combined['Air_Visitors'] + df_combined['Cruise_Visitors']
            return df_combined

    # Fallback: air visitors only
    # BUG FIX: df_air is a DataFrame, not a dict — use bracket access or .get on columns
    df_air['Cruise_Visitors'] = 0
    df_air['Total_Bonaire_Visitors'] = df_air['Air_Visitors']
    return df_air


# -----------------------------------------------------------------------------
# 3. CDSE API FUNCTIONS: CATALOGUE (DATES) & SENTINEL HUB (PROCESSING)
# -----------------------------------------------------------------------------
@st.cache_data(ttl=3600)
def get_available_dates(sensor_collection, year, month):
    """
    Queries CDSE OData API to find exactly which days have satellite overpasses.

    BUG FIX: The CDSE catalogue collection name is 'SENTINEL-2' (uppercase), NOT
    'sentinel-2-l2a' or 'sentinel-2-l1c' (those are Sentinel Hub processing API IDs).
    The L1C/L2A product type is a separate Attributes filter.
    """
    last_day = calendar.monthrange(year, month)[1]
    start_date = f"{year}-{month:02d}-01T00:00:00.000Z"
    end_date   = f"{year}-{month:02d}-{last_day}T23:59:59.000Z"

    # Map Sentinel Hub collection IDs -> CDSE catalogue names + product type codes
    COLLECTION_MAP = {
        "sentinel-2-l2a": ("SENTINEL-2", "S2MSI2A"),
        "sentinel-2-l1c": ("SENTINEL-2", "S2MSI1C"),
    }
    cdse_collection, product_type = COLLECTION_MAP.get(
        sensor_collection, ("SENTINEL-2", "S2MSI2A")
    )

    bbox = "POLYGON((-68.3 12.0, -68.1 12.0, -68.1 12.3, -68.3 12.3, -68.3 12.0))"

    # Correct CDSE OData filter with productType attribute
    url = (
        f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
        f"?$filter=Collection/Name eq '{cdse_collection}'"
        f" and Attributes/OData.CSC.StringAttribute/any("
        f"att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq '{product_type}')"
        f" and OData.CSC.Intersects(area=geography'SRID=4326;{bbox}')"
        f" and ContentDate/Start ge {start_date}"
        f" and ContentDate/Start le {end_date}"
        f"&$select=ContentDate&$top=100"
    )

    try:
        resp = requests.get(url, timeout=20)
        if resp.status_code != 200:
            st.error(f"CDSE Catalogue HTTP {resp.status_code}: {resp.text[:300]}")
            return []
        data = resp.json()
        if 'value' in data:
            dates = [item['ContentDate']['Start'][:10] for item in data['value']]
            return sorted(list(set(dates)))
        else:
            st.error(f"Unexpected CDSE response structure: {str(data)[:300]}")
    except Exception as e:
        st.error(f"Catalogue request error: {e}")
    return []


def get_sh_token(client_id, client_secret):
    """Obtains a Sentinel Hub OAuth2 token."""
    url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    payload = {
        "grant_type": "client_credentials",
        "client_id": client_id,
        "client_secret": client_secret,
    }
    try:
        response = requests.post(url, data=payload, timeout=15)
        if response.status_code == 200:
            return response.json().get("access_token")
        else:
            st.error(f"Token request failed ({response.status_code}): {response.text}")
    except Exception as e:
        st.error(f"Token request error: {e}")
    return None


def build_evalscript(fai_threshold, mci_threshold):
    """
    Builds a Sentinel Hub evalscript that:
      - Computes FAI (Floating Algae Index) and MCI (Maximum Chlorophyll Index)
      - Flags pixels exceeding both thresholds as sargassum
      - Returns a 3-band image: [sargassum_mask, FAI, MCI]
    """
    return f"""
//VERSION=3
function setup() {{
  return {{
    input: [{{
      bands: ["B04", "B06", "B08", "B8A", "B11"],
      units: "REFLECTANCE"
    }}],
    output: {{
      bands: 3,
      sampleType: "FLOAT32"
    }}
  }};
}}

function evaluatePixel(sample) {{
  // FAI: Floating Algae Index
  // Reference: Hu (2009), doi:10.1029/2009JC005370
  var fai = sample.B08 - (sample.B04 + (sample.B11 - sample.B04) * ((832.8 - 664.6) / (1613.7 - 664.6)));

  // MCI: Maximum Chlorophyll Index
  // Reference: Gower et al. (2005)
  var mci = sample.B8A - sample.B06 - (sample.B11 - sample.B06) * ((864.8 - 740.5) / (1613.7 - 740.5));

  // Sargassum mask: both indices must exceed thresholds
  var sargassum = (fai > {fai_threshold} && mci > {mci_threshold}) ? 1.0 : 0.0;

  return [sargassum, fai, mci];
}}
"""


def process_sargassum_image(token, date_val, collection, evalscript):
    """
    Calls the Sentinel Hub Processing API for Bonaire's bounding box.
    Returns a PIL Image (RGB) or None on failure.
    BUG FIX: completed the truncated payload — 'to' field and all closing braces were missing.
    """
    url = "https://sh.dataspace.copernicus.eu/api/v1/process"
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json"
    }
    payload = {
        "input": {
            "bounds": {
                "bbox": [-68.3, 12.0, -68.1, 12.3],
                "properties": {"crs": "http://www.opengis.net/def/crs/EPSG/0/4326"}
            },
            "data": [{
                "type": collection,
                "dataFilter": {
                    "timeRange": {
                        "from": f"{date_val}T00:00:00Z",
                        "to":   f"{date_val}T23:59:59Z"   # BUG FIX: was missing entirely
                    }
                }
            }]
        },
        "output": {
            "width": 512,
            "height": 512,
            "responses": [{
                "identifier": "default",
                "format": {"type": "image/tiff"}
            }]
        },
        "evalscript": evalscript
    }
    try:
        response = requests.post(url, headers=headers, json=payload, timeout=60)
        if response.status_code == 200:
            return Image.open(io.BytesIO(response.content))
        else:
            st.error(f"Processing API error ({response.status_code}): {response.text[:300]}")
    except Exception as e:
        st.error(f"Processing request failed: {e}")
    return None


def estimate_sargassum_area_km2(img):
    """
    Given a PIL image where band 0 is the sargassum mask (1.0 = sargassum),
    estimates the covered area in km².
    Bonaire bbox: lon [-68.3, -68.1] x lat [12.0, 12.3] ≈ 22.2 km × 33.3 km = ~740 km²
    """
    arr = np.array(img)
    # Handle both single-band and multi-band TIFF
    if arr.ndim == 3:
        mask = arr[:, :, 0]
    else:
        mask = arr
    total_pixels = mask.size
    sarg_pixels  = np.sum(mask > 0.5)
    bbox_area_km2 = 22.2 * 33.3  # approximate
    return (sarg_pixels / total_pixels) * bbox_area_km2


# -----------------------------------------------------------------------------
# 4. REGRESSION
# -----------------------------------------------------------------------------
def run_regression(df_merged, y_col, x_col, lag_months=0):
    """
    OLS regression of y_col ~ x_col with optional lag on x.
    Returns the statsmodels RegressionResults object.
    """
    df = df_merged[[y_col, x_col]].dropna().copy()
    if lag_months > 0:
        df[x_col] = df[x_col].shift(lag_months)
        df = df.dropna()
    if len(df) < 5:
        return None
    X = sm.add_constant(df[x_col])
    model = sm.OLS(df[y_col], X).fit()
    return model, df


# -----------------------------------------------------------------------------
# 5. STREAMLIT UI
# -----------------------------------------------------------------------------

# ── Sidebar: credentials & settings ──────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Settings")

    st.subheader("CDSE Credentials")
    client_id     = st.text_input("Client ID",     type="password")
    client_secret = st.text_input("Client Secret", type="password")

    st.subheader("Sensor")
    collection = st.selectbox(
        "Sentinel Collection",
        ["sentinel-2-l2a", "sentinel-2-l1c"],
        help="L2A (surface reflectance) is recommended for water applications."
    )

    st.subheader("Detection Thresholds")
    fai_threshold = st.slider("FAI Threshold", min_value=-0.05, max_value=0.10, value=0.01, step=0.005,
                               help="Floating Algae Index. Increase to reduce false positives.")
    mci_threshold = st.slider("MCI Threshold", min_value=-0.01, max_value=0.05, value=0.005, step=0.001,
                               help="Maximum Chlorophyll Index.")

    st.subheader("Date Selection")
    year  = st.selectbox("Year",  list(range(2017, datetime.date.today().year + 1))[::-1])
    month = st.selectbox("Month", list(range(1, 13)),
                          format_func=lambda m: calendar.month_name[m])

    load_dates_btn = st.button("🔍 Find Available Satellite Dates")

# ── Main tabs ─────────────────────────────────────────────────────────────────
tab_sat, tab_tourism, tab_regression = st.tabs(
    ["🛰️ Satellite Processing", "📊 Tourism Data", "📈 Regression Analysis"]
)

# ── Tab 1: Satellite ──────────────────────────────────────────────────────────
with tab_sat:
    st.header("Sargassum Detection — Bonaire")

    available_dates = []
    if load_dates_btn:
        if not client_id or not client_secret:
            st.warning("Please enter your CDSE credentials in the sidebar.")
        else:
            with st.spinner("Querying CDSE catalogue..."):
                available_dates = get_available_dates(collection, year, month)
            if available_dates:
                st.success(f"Found {len(available_dates)} overpass(es) for {calendar.month_name[month]} {year}.")
            else:
                st.warning("No satellite data found for this month/area. Try a different month.")

    if available_dates:
        selected_date = st.selectbox(
            "Select acquisition date",
            available_dates,
            help="Only dates with confirmed satellite coverage over Bonaire are shown."
        )
        process_btn = st.button("🚀 Process Image & Detect Sargassum")

        if process_btn:
            with st.spinner("Authenticating with Sentinel Hub..."):
                token = get_sh_token(client_id, client_secret)

            if token:
                evalscript = build_evalscript(fai_threshold, mci_threshold)
                with st.spinner(f"Processing {selected_date}... (this may take ~30 s)"):
                    img = process_sargassum_image(token, selected_date, collection, evalscript)

                if img is not None:
                    arr = np.array(img)
                    # Display the sargassum mask
                    col1, col2 = st.columns(2)
                    with col1:
                        st.subheader("Sargassum Mask")
                        if arr.ndim == 3:
                            mask_display = arr[:, :, 0]
                        else:
                            mask_display = arr
                        fig, ax = plt.subplots(figsize=(5, 5))
                        ax.imshow(mask_display, cmap='YlGn', vmin=0, vmax=1)
                        ax.set_title(f"Sargassum Mask — {selected_date}")
                        ax.axis('off')
                        st.pyplot(fig)

                    with col2:
                        area_km2 = estimate_sargassum_area_km2(img)
                        st.metric("Estimated Sargassum Coverage", f"{area_km2:.2f} km²")
                        st.info(
                            f"FAI threshold: {fai_threshold}  |  MCI threshold: {mci_threshold}\n\n"
                            "Adjust thresholds in the sidebar to tune detection sensitivity."
                        )

                    # Store result in session state for regression tab
                    if 'sargassum_records' not in st.session_state:
                        st.session_state['sargassum_records'] = []
                    st.session_state['sargassum_records'].append({
                        'Date': pd.to_datetime(selected_date),
                        'Sargassum_km2': area_km2
                    })
                    st.success("Result saved. Switch to the Regression tab to analyse it against tourism data.")
            else:
                st.error("Authentication failed. Check your credentials.")
    else:
        st.info("Use the sidebar to select a year/month and click **Find Available Satellite Dates**.")

# ── Tab 2: Tourism Data ───────────────────────────────────────────────────────
with tab_tourism:
    st.header("Tourism Data — Bonaire (CBS / Statistics Netherlands)")

    load_tourism_btn = st.button("📥 Load Tourism Data")

    if load_tourism_btn or 'tourism_data' in st.session_state:
        if load_tourism_btn:
            with st.spinner("Fetching data from CBS Open Data API..."):
                df_tourism = load_tourism_data()
            st.session_state['tourism_data'] = df_tourism
        else:
            df_tourism = st.session_state['tourism_data']

        if df_tourism.empty:
            st.error("No tourism data could be loaded. Check CBS table availability.")
        else:
            st.success(f"Loaded {len(df_tourism)} monthly records "
                       f"({df_tourism.index.min().strftime('%Y-%m')} – {df_tourism.index.max().strftime('%Y-%m')}).")

            st.subheader("Monthly Visitor Counts")
            fig, ax = plt.subplots(figsize=(12, 4))
            if 'Air_Visitors' in df_tourism.columns:
                ax.bar(df_tourism.index, df_tourism['Air_Visitors'], label='Air', alpha=0.8)
            if 'Cruise_Visitors' in df_tourism.columns:
                ax.bar(df_tourism.index, df_tourism['Cruise_Visitors'],
                       bottom=df_tourism.get('Air_Visitors', 0), label='Cruise', alpha=0.8)
            ax.set_xlabel("Date")
            ax.set_ylabel("Visitors")
            ax.legend()
            ax.set_title("Bonaire — Monthly Inbound Visitors")
            st.pyplot(fig)

            with st.expander("View raw data"):
                st.dataframe(df_tourism)
    else:
        st.info("Click **Load Tourism Data** to fetch from Statistics Netherlands (CBS).")

# ── Tab 3: Regression ─────────────────────────────────────────────────────────
with tab_regression:
    st.header("Regression: Sargassum Influx vs. Tourism")

    # Build sargassum monthly series from session state
    if 'sargassum_records' in st.session_state and st.session_state['sargassum_records']:
        df_sarg = pd.DataFrame(st.session_state['sargassum_records'])
        df_sarg['Month'] = df_sarg['Date'].dt.to_period('M').dt.to_timestamp()
        df_sarg = df_sarg.groupby('Month')['Sargassum_km2'].mean().to_frame()
        df_sarg.index.name = 'Date'
        st.info(f"Sargassum records in session: {len(df_sarg)} month(s). "
                "Process more dates in the Satellite tab to build a longer series.")
    else:
        st.warning("No sargassum data yet. Process satellite images in the **Satellite Processing** tab first.")
        df_sarg = pd.DataFrame()

    if 'tourism_data' in st.session_state and not st.session_state['tourism_data'].empty:
        df_tourism = st.session_state['tourism_data']
    else:
        st.info("Load tourism data in the **Tourism Data** tab first.")
        df_tourism = pd.DataFrame()

    if not df_sarg.empty and not df_tourism.empty:
        # Merge on monthly date index
        df_merged = df_tourism.join(df_sarg, how='inner')

        if df_merged.empty:
            st.warning("No overlapping dates between sargassum records and tourism data.")
        else:
            st.subheader("Merged Dataset")
            st.dataframe(df_merged)

            y_col = st.selectbox("Dependent variable (Y)",
                                  [c for c in df_merged.columns if 'Visitor' in c or 'visitor' in c],
                                  help="Tourism metric to predict.")
            lag = st.slider("Lag (months)", 0, 6, 0,
                            help="Shift sargassum data forward by N months before regressing.")

            run_reg_btn = st.button("▶ Run OLS Regression")
            if run_reg_btn:
                result = run_regression(df_merged, y_col, 'Sargassum_km2', lag_months=lag)
                if result is None:
                    st.error("Not enough data points after applying lag. Process more satellite dates.")
                else:
                    model, df_reg = result
                    st.subheader("Regression Results")
                    st.text(model.summary().as_text())

                    # Scatter + regression line
                    fig2, ax2 = plt.subplots(figsize=(8, 5))
                    ax2.scatter(df_reg['Sargassum_km2'], df_reg[y_col], alpha=0.7, label='Observations')
                    x_line = np.linspace(df_reg['Sargassum_km2'].min(), df_reg['Sargassum_km2'].max(), 100)
                    y_line = model.params['const'] + model.params['Sargassum_km2'] * x_line
                    ax2.plot(x_line, y_line, color='red', label='OLS fit')
                    ax2.set_xlabel(f"Sargassum Coverage (km²){f'  [lag {lag}m]' if lag else ''}")
                    ax2.set_ylabel(y_col)
                    ax2.set_title(f"OLS: {y_col} ~ Sargassum_km2  |  R²={model.rsquared:.3f}  p={model.pvalues['Sargassum_km2']:.4f}")
                    ax2.legend()
                    st.pyplot(fig2)

                    # Aruba comparison note
                    st.info(
                        "💡 **Tip (Aruba control):** Your proposal suggests using Aruba as a control island. "
                        "Load Aruba's CBS data (same table IDs, filter for 'Aruba') and compare the sargassum "
                        "coefficient between islands to infer causality."
                    )
    else:
        st.info("Complete both the Satellite Processing and Tourism Data tabs to enable regression.")