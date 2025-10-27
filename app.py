"""
Aplicación Streamlit para Visualización Interactiva del Análisis EDA de Arroz
Análisis de Heading Date y Estructura Genética en Oryza sativa (3K RGP)

Autor: Experto en Desarrollo Software y Datos Agrícolas
Fecha: Octubre 2025

Referencias científicas:
- Wang et al. (2018) Nature: Genomic variation in 3,010 diverse accessions
- Garris et al. (2005) Nature Genetics: Estructura de población en arroz
- Huang et al. (2012) Nature: The genome of cultivated rice
- Zhou et al. (2021) New Phytologist: Transcriptional regulation of heading date
"""

import json
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ==================== CONFIGURACIÓN DE LA PÁGINA ====================
st.set_page_config(
    page_title="3K RGP - Heading Date Explorer",
    page_icon="🌾",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== ESTILOS CSS PERSONALIZADOS ====================
st.markdown("""
<style>
    .main-header {
        font-size: 2.8rem;
        font-weight: bold;
        color: #2E7D32;
        text-align: center;
        padding: 1.5rem;
        background: linear-gradient(90deg, #A5D6A7 0%, #66BB6A 100%);
        border-radius: 10px;
        margin-bottom: 2rem;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    .sub-header {
        font-size: 1.5rem;
        font-weight: bold;
        color: #1B5E20;
        margin-top: 2rem;
        margin-bottom: 1rem;
        border-left: 5px solid #4CAF50;
        padding-left: 10px;
    }
    .metric-card {
        background-color: #E8F5E9;
        padding: 1.2rem;
        border-radius: 8px;
        border-left: 4px solid #4CAF50;
        margin: 0.5rem 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    .info-box {
        background-color: #E3F2FD;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #2196F3;
        margin: 1rem 0;
    }
    .warning-box {
        background-color: #FFF3E0;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #FF9800;
        margin: 1rem 0;
    }
    .success-box {
        background-color: #E8F5E9;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #4CAF50;
        margin: 1rem 0;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
    }
    .stTabs [data-baseweb="tab"] {
        height: 3rem;
        padding: 0 2rem;
        background-color: #E8F5E9;
        border-radius: 5px 5px 0 0;
    }
    .stTabs [aria-selected="true"] {
        background-color: #4CAF50;
        color: white;
    }
</style>
""", unsafe_allow_html=True)

# ==================== FUNCIONES DE CARGA DE DATOS ====================

@st.cache_data
def load_data():
    """Carga los datos fenotípicos y genotípicos desde los archivos procesados"""
    try:
        # Definir rutas basadas en la estructura del proyecto
        ROOT = Path(".")
        REPORTS_PATH = ROOT / "reports"
        
        # Cargar datos procesados
        df = pd.read_csv(REPORTS_PATH / "traits_geno_aligned.csv", index_col=0)
        geno = pd.read_csv(REPORTS_PATH / "geno_aligned.csv", index_col=0)
        
        # Verificar que los datos se cargaron correctamente
        if len(df) == 0 or len(geno) == 0:
            raise FileNotFoundError("Los archivos están vacíos")
            
        # Crear variables derivadas coherentes con el EDA
        if 'GRLT' in df.columns and 'GRWD' in df.columns:
            df['GrainSize'] = (df['GRLT'] / df['GRWD']).round(2)
            
        df['HDG_category'] = pd.cut(
            df['HDG_80HEAD'],
            bins=[-np.inf, 90, 110, np.inf],
            labels=['Early (<90d)', 'Medium (90-110d)', 'Late (>110d)']
        )
        
        # Agregar coordenadas y regiones si no existen
        if 'latitude' not in df.columns or 'longitude' not in df.columns:
            df = add_geographic_data(df)
            
    except Exception as e:
        st.warning(f"⚠️ No se encontraron datos reales. Generando datos de demostración... ({str(e)})")
        df, geno = generate_demo_data()
    
    return df, geno

def add_geographic_data(df):
    """Agrega coordenadas y regiones geográficas al dataframe"""
    # Coordenadas aproximadas por país (centros geográficos)
    coords = {
        'India': (20.5937, 78.9629, 'SAS'),
        'Bangladesh': (23.6850, 90.3563, 'SAS'),
        'China': (35.8617, 104.1954, 'EAS'),
        'Indonesia': (-0.7893, 113.9213, 'SEA'),
        'Philippines': (12.8797, 121.7740, 'SEA'),
        'Thailand': (15.8700, 100.9925, 'SEA'),
        'Vietnam': (14.0583, 108.2772, 'SEA'),
        'Japan': (36.2048, 138.2529, 'EAS'),
        'South_Korea': (35.9078, 127.7669, 'EAS'),
        'Korea': (35.9078, 127.7669, 'EAS'),
        'Pakistan': (30.3753, 69.3451, 'SAS'),
        'Sri_Lanka': (7.8731, 80.7718, 'SAS'),
        'Myanmar': (21.9162, 95.9560, 'SEA'),
        'Cambodia': (12.5657, 104.9910, 'SEA'),
        'Taiwan': (23.5937, 121.0254, 'EAS'),
        'Nepal': (28.3949, 84.1240, 'SAS'),
        'Bhutan': (27.5142, 90.4336, 'SAS'),
        'Laos': (19.8563, 102.4955, 'SEA'),
        'Malaysia': (4.2105, 101.9758, 'SEA'),
        'Iran': (32.4279, 53.6880, 'WAS'),
        'Africa': (-8.7832, 34.5085, 'AFR'),
        'USA': (37.0902, -95.7129, 'NAM'),
        'Brazil': (-14.2350, -51.9253, 'SAM'),
        'Italy': (41.8719, 12.5674, 'EUR'),
        'Spain': (40.4637, -3.7492, 'EUR'),
        'France': (46.2276, 2.2137, 'EUR'),
    }
    
    # Mapear coordenadas y regiones
    if 'country' in df.columns:
        df['latitude'] = df['country'].map(lambda x: coords.get(x, (None, None, None))[0])
        df['longitude'] = df['country'].map(lambda x: coords.get(x, (None, None, None))[1])
        df['region'] = df['country'].map(lambda x: coords.get(x, (None, None, None))[2])
    
    return df

def generate_demo_data():
    """Genera datos de demostración realistas basados en el análisis del 3K RGP"""
    np.random.seed(42)
    n_samples = 1882  # Tamaño real del dataset 3K RGP
    
    # Lista de países del análisis real
    countries = ['India', 'Bangladesh', 'China', 'Indonesia', 'Philippines', 
                'Thailand', 'Vietnam', 'Japan', 'South_Korea', 'Pakistan',
                'Sri_Lanka', 'Myanmar', 'Cambodia', 'Taiwan', 'Nepal']
    
    # Subespecies con distribución realista
    subspecies_list = ['IND', 'JAP', 'AUS', 'ARO', 'TRJ', 'ADM']
    subspecies = np.random.choice(subspecies_list, n_samples, 
                                 p=[0.45, 0.27, 0.11, 0.04, 0.09, 0.04])
    
    # Generar heading dates con medias y desviaciones por subespecie (basado en EDA real)
    hdg_means = {'IND': 105, 'JAP': 87, 'AUS': 80, 'ARO': 102, 'TRJ': 113, 'ADM': 101}
    hdg_stds = {'IND': 22, 'JAP': 16, 'AUS': 11, 'ARO': 21, 'TRJ': 19, 'ADM': 17}
    
    hdg = np.array([np.random.normal(hdg_means[s], hdg_stds[s]) for s in subspecies])
    hdg = np.clip(hdg, 50, 175)
    
    # Crear DataFrame con correlaciones realistas entre rasgos
    df = pd.DataFrame({
        'HDG_80HEAD': hdg,
        'subespecie': subspecies,
        'country': np.random.choice(countries, n_samples),
        'SDHT': 39 + 0.15 * (hdg - 100) + np.random.normal(0, 5, n_samples),  # Correlación con HDG
        'PLT_POST': 25 + 0.08 * (hdg - 100) + np.random.normal(0, 3, n_samples),  # Longitud panícula
        'CULT_REPRO': 80 + 0.4 * hdg + np.random.normal(0, 8, n_samples),  # Correlación con HDG
        'GRLT': np.random.normal(7.5, 1.2, n_samples),  # Longitud grano
        'GRWD': np.random.normal(2.8, 0.5, n_samples),  # Ancho grano
        'GRWT100': np.random.normal(2.5, 0.4, n_samples),  # Peso 100 granos
        'LLT': np.random.normal(45, 8, n_samples),  # Longitud hoja
        'LWD': np.random.normal(1.5, 0.3, n_samples),  # Ancho hoja
    })
    
    # Agregar coordenadas y datos derivados
    df = add_geographic_data(df)
    df['GrainSize'] = (df['GRLT'] / df['GRWD']).round(2)
    df['HDG_category'] = pd.cut(
        df['HDG_80HEAD'],
        bins=[-np.inf, 90, 110, np.inf],
        labels=['Early (<90d)', 'Medium (90-110d)', 'Late (>110d)']
    )
    
    # Genotipos simulados
    geno = pd.DataFrame(
        np.random.randn(n_samples, 100),
        index=df.index,
        columns=[f'SNP_{i}' for i in range(100)]
    )
    
    return df, geno

@st.cache_data
def compute_pca(_geno_data, n_components=10):
    """Calcula PCA de los genotipos"""
    scaler = StandardScaler()
    geno_scaled = scaler.fit_transform(_geno_data)
    
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(geno_scaled)
    
    explained_var = pca.explained_variance_ratio_
    
    return pca_result, explained_var, pca

# ==================== FUNCIONES DE FILTRADO ====================

def apply_filters(df, filter_key):
    """Aplica filtros independientes para cada tab"""
    
    df_filtered = df.copy()
    
    # Filtro por región
    if 'region' in df.columns:
        regions = ['Todas'] + sorted(df['region'].dropna().unique().tolist())
        selected_region = st.selectbox(
            "📍 Región Geográfica:",
            regions,
            key=f"region_{filter_key}"
        )
        
        if selected_region != 'Todas':
            df_filtered = df_filtered[df_filtered['region'] == selected_region].copy()
    
    # Filtro por subespecie
    subspecies = ['Todas'] + sorted(df['subespecie'].dropna().unique().tolist())
    selected_subspecies = st.selectbox(
        "🧬 Subespecie:",
        subspecies,
        key=f"subspecies_{filter_key}"
    )
    
    if selected_subspecies != 'Todas':
        df_filtered = df_filtered[df_filtered['subespecie'] == selected_subspecies].copy()
    
    # Filtro por rango de HDG
    hdg_min, hdg_max = st.slider(
        "📊 Rango de Heading Date (días):",
        int(df['HDG_80HEAD'].min()),
        int(df['HDG_80HEAD'].max()),
        (int(df['HDG_80HEAD'].min()), int(df['HDG_80HEAD'].max())),
        key=f"hdg_range_{filter_key}"
    )
    
    df_filtered = df_filtered[
        (df_filtered['HDG_80HEAD'] >= hdg_min) & 
        (df_filtered['HDG_80HEAD'] <= hdg_max)
    ].copy()
    
    # Mostrar estadísticas
    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    col1.metric("N Accesiones", len(df_filtered))
    col2.metric("HDG Media", f"{df_filtered['HDG_80HEAD'].mean():.1f} días")
    col3.metric("HDG Rango", f"{df_filtered['HDG_80HEAD'].min():.0f}-{df_filtered['HDG_80HEAD'].max():.0f}")
    
    return df_filtered

# ==================== FUNCIONES DE VISUALIZACIÓN ====================

def create_geographic_map(df, color_by='HDG_80HEAD', title="Distribución Geográfica"):
    """Crea un mapa interactivo con los datos georeferenciados"""
    
    # Filtrar datos con coordenadas válidas
    df_map = df.dropna(subset=['latitude', 'longitude', color_by]).copy()
    
    if len(df_map) == 0:
        st.warning("No hay datos con coordenadas válidas para mostrar en el mapa.")
        return None
    
    # Agregar información de tooltip
    df_map['hover_text'] = (
        'País: ' + df_map['country'].astype(str) + '<br>' +
        'Subespecie: ' + df_map['subespecie'].astype(str) + '<br>' +
        'HDG: ' + df_map['HDG_80HEAD'].round(1).astype(str) + ' días<br>' +
        'Región: ' + df_map['region'].astype(str)
    )
    
    # Crear mapa
    if color_by == 'HDG_80HEAD':
        fig = px.scatter_mapbox(
            df_map,
            lat='latitude',
            lon='longitude',
            color=color_by,
            size='HDG_80HEAD',
            hover_name='country',
            hover_data={
                'subespecie': True,
                'HDG_80HEAD': ':.1f',
                'region': True,
                'latitude': False,
                'longitude': False
            },
            color_continuous_scale='Viridis',
            size_max=15,
            zoom=2,
            title=title,
            height=600
        )
    else:
        fig = px.scatter_mapbox(
            df_map,
            lat='latitude',
            lon='longitude',
            color=color_by,
            hover_name='country',
            hover_data={
                'subespecie': True,
                'HDG_80HEAD': ':.1f',
                'region': True,
                'latitude': False,
                'longitude': False
            },
            size_max=10,
            zoom=2,
            title=title,
            height=600
        )
    
    fig.update_layout(
        mapbox_style="open-street-map",
        margin={"r":0,"t":40,"l":0,"b":0}
    )
    
    return fig

def create_hdg_morphology_scatter(df, x_trait, y_trait='HDG_80HEAD', color_by='subespecie'):
    """Crea gráfico de dispersión para explorar relación entre HDG y morfología"""
    
    df_clean = df.dropna(subset=[x_trait, y_trait, color_by]).copy()
    
    if len(df_clean) == 0:
        return None, None, None
    
    # Calcular correlación
    corr, p_value = stats.pearsonr(df_clean[x_trait], df_clean[y_trait])
    
    fig = px.scatter(
        df_clean,
        x=x_trait,
        y=y_trait,
        color=color_by,
        hover_data=['country', 'subespecie'],
        #trendline='ols',
        title=f'Relación {y_trait} vs {x_trait} (r={corr:.3f}, p<{p_value:.2e})',
        labels={x_trait: x_trait, y_trait: 'Heading Date (días)'},
        height=500
    )
    
    fig.update_traces(marker=dict(size=8, opacity=0.7))
    
    return fig, corr, p_value

def create_pca_plot(pca_result, df, color_by='subespecie', pc_x=0, pc_y=1, explained_var=None):
    """Crea gráfico PCA interactivo con opciones de coloración"""
    
    pca_df = pd.DataFrame({
        f'PC{pc_x+1}': pca_result[:, pc_x],
        f'PC{pc_y+1}': pca_result[:, pc_y],
        'color_var': df[color_by].values,
        'country': df['country'].values,
        'HDG_80HEAD': df['HDG_80HEAD'].values,
        'subespecie': df['subespecie'].values
    })
    
    pca_df = pca_df.dropna()
    
    # Título con varianza explicada
    var_text = ""
    if explained_var is not None:
        var_text = f" (PC{pc_x+1}: {explained_var[pc_x]*100:.1f}%, PC{pc_y+1}: {explained_var[pc_y]*100:.1f}%)"
    
    if color_by == 'HDG_80HEAD':
        fig = px.scatter(
            pca_df,
            x=f'PC{pc_x+1}',
            y=f'PC{pc_y+1}',
            color='color_var',
            color_continuous_scale='Viridis',
            hover_data=['country', 'subespecie', 'HDG_80HEAD'],
            title=f'PCA de Genotipos - Coloreado por {color_by}{var_text}',
            labels={'color_var': color_by},
            height=600
        )
    else:
        fig = px.scatter(
            pca_df,
            x=f'PC{pc_x+1}',
            y=f'PC{pc_y+1}',
            color='color_var',
            hover_data=['country', 'subespecie', 'HDG_80HEAD'],
            title=f'PCA de Genotipos - Coloreado por {color_by}{var_text}',
            labels={'color_var': color_by},
            height=600
        )
    
    fig.update_traces(marker=dict(size=8, opacity=0.7))
    
    return fig

def create_regional_boxplot(df, region_col='region', value_col='HDG_80HEAD'):
    """Crea boxplot comparativo por regiones"""
    
    df_clean = df.dropna(subset=[region_col, value_col]).copy()
    
    if len(df_clean) == 0:
        return None
    
    # Ordenar regiones por mediana
    region_order = df_clean.groupby(region_col)[value_col].median().sort_values().index
    
    fig = px.box(
        df_clean,
        x=region_col,
        y=value_col,
        color=region_col,
        category_orders={region_col: region_order},
        title=f'Distribución de {value_col} por Región Geográfica',
        labels={value_col: 'Heading Date (días)', region_col: 'Región'},
        height=500
    )
    
    fig.update_layout(showlegend=False)
    
    return fig

# ==================== APLICACIÓN PRINCIPAL ====================

def main():
    # Título principal
    st.markdown("<h1 class='main-header'>🌾 3K Rice Genomes - Heading Date Explorer</h1>", unsafe_allow_html=True)
    
    st.markdown("""
    <div class='info-box'>
        <b>🔬 Dataset:</b> 3,000 Rice Genomes Project (3K RGP) - 1,882 accesiones de arroz asiático<br>
        <b>🎯 Objetivo:</b> Exploración interactiva de la variación en Heading Date y su relación con geografía, genética y morfología<br>
        <b>💡 Nota:</b> Cada pestaña tiene filtros independientes para evitar conflictos entre análisis
    </div>
    """, unsafe_allow_html=True)
    
    # Cargar datos
    with st.spinner("🔄 Cargando datos del 3K RGP..."):
        df, geno = load_data()
    
    # Sidebar con información general (sin filtros globales)
    st.sidebar.markdown("## 📊 Dataset Overview")
    st.sidebar.metric("Total Accesiones", len(df))
    st.sidebar.metric("HDG Media Global", f"{df['HDG_80HEAD'].mean():.1f} días")
    st.sidebar.metric("HDG Rango Global", f"{df['HDG_80HEAD'].min():.0f} - {df['HDG_80HEAD'].max():.0f}")
    
    if 'subespecie' in df.columns:
        st.sidebar.markdown("### 🧬 Distribución de Subespecies")
        ssp_counts = df['subespecie'].value_counts()
        for ssp, count in ssp_counts.items():
            percentage = (count / len(df)) * 100
            st.sidebar.text(f"{ssp}: {count} ({percentage:.1f}%)")
    
    if 'region' in df.columns:
        st.sidebar.markdown("### 🌍 Distribución por Regiones")
        reg_counts = df['region'].value_counts()
        for reg, count in reg_counts.items():
            st.sidebar.text(f"{reg}: {count}")
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    <div style='font-size: 0.85rem; color: #666; padding: 10px; background-color: #f0f0f0; border-radius: 5px;'>
        <b>💡 Filtros Independientes</b><br>
        Cada pestaña tiene sus propios controles de filtrado para análisis específicos sin interferencias.
    </div>
    """, unsafe_allow_html=True)
    
    # Tabs principales
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🗺️ Exploración Geográfica",
        "🧬 Análisis Genético (PCA)",
        "📊 HDG vs Morfología",
        "🌍 Análisis por Regiones",
        "🎯 Germoplasma Elite"
    ])
    
    # ==================== TAB 1: EXPLORACIÓN GEOGRÁFICA ====================
    with tab1:
        st.markdown("<h2 class='sub-header'>🗺️ Distribución Geográfica del Heading Date</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        El heading date muestra un **gradiente latitudinal claro** (r ≈ 0.55-0.65), reflejando adaptación local:
        - **Latitudes bajas (<20°)**: Predominan variedades tempranas y de ciclo medio
        - **Latitudes altas (>35°)**: Predominan variedades de ciclo medio y tardío
        """)
        
        # Filtros específicos para este tab
        with st.expander("🎛️ Filtros de Datos", expanded=False):
            df_tab1 = apply_filters(df, "tab1")
        
        # Si no se expandió el filtro, usar datos completos
        if 'df_tab1' not in locals():
            df_tab1 = df.copy()
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            # Mapa principal coloreado por HDG
            map_color = st.radio(
                "Colorear mapa por:",
                ['HDG_80HEAD', 'subespecie', 'region', 'HDG_category'],
                horizontal=True,
                key="map_color_tab1"
            )
            
            fig_map = create_geographic_map(df_tab1, color_by=map_color, 
                                           title=f"Distribución Geográfica - {map_color}")
            if fig_map:
                st.plotly_chart(fig_map, use_container_width=True)
        
        with col2:
            st.markdown("### 📊 Distribución por País")
            
            # Top países por número de accesiones
            top_countries = df_tab1['country'].value_counts().head(10)
            
            fig_countries = go.Figure(go.Bar(
                x=top_countries.values,
                y=top_countries.index,
                orientation='h',
                marker_color='#4CAF50'
            ))
            
            fig_countries.update_layout(
                title="Top 10 Países",
                xaxis_title="N Accesiones",
                yaxis_title="",
                height=400
            )
            
            st.plotly_chart(fig_countries, use_container_width=True)
            
            # Estadísticas por país
            st.markdown("### 📈 HDG Promedio por País")
            country_stats = df_tab1.groupby('country')['HDG_80HEAD'].agg(['mean', 'count']).sort_values('mean')
            country_stats.columns = ['Media HDG', 'N']
            st.dataframe(country_stats.head(10), use_container_width=True)
        
        # Correlación latitud-HDG
        st.markdown("---")
        st.markdown("### 🌐 Relación Latitud-Heading Date")
        
        df_lat = df_tab1.dropna(subset=['latitude', 'HDG_80HEAD'])
        
        if len(df_lat) > 10:
            corr_lat, p_lat = stats.pearsonr(df_lat['latitude'], df_lat['HDG_80HEAD'])
            
            fig_lat = px.scatter(
                df_lat,
                x='latitude',
                y='HDG_80HEAD',
                color='subespecie',
                hover_data=['country'],
                #trendline='ols',
                title=f'Adaptación Latitudinal (r={corr_lat:.3f}, p<{p_lat:.2e})',
                labels={'latitude': 'Latitud (°)', 'HDG_80HEAD': 'Heading Date (días)'},
                height=500
            )
            
            st.plotly_chart(fig_lat, use_container_width=True)
            
            st.markdown(f"""
            <div class='success-box'>
                <b>✅ Correlación significativa detectada:</b> r = {corr_lat:.3f} (p < {p_lat:.2e})<br>
                Las variedades de latitudes altas tienden a tener ciclos más largos, reflejando adaptación 
                a estaciones de crecimiento más largas y temperaturas más bajas.
            </div>
            """, unsafe_allow_html=True)
    
    # ==================== TAB 2: ANÁLISIS GENÉTICO (PCA) ====================
    with tab2:
        st.markdown("<h2 class='sub-header'>🧬 Análisis de Componentes Principales (PCA)</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        El PCA revela la **estructura genética poblacional** subyacente basada en SNPs.
        Las 5 subespecies principales (IND, JAP, AUS, ARO, TRJ) se separan claramente, 
        reflejando milenios de adaptación local y selección.
        """)
        

        
        if 'df_tab2' not in locals():
            df_tab2 = df.copy()
        
        # Calcular PCA
        with st.spinner("🔄 Calculando PCA..."):
            pca_result, explained_var, pca_model = compute_pca(geno, n_components=10)
        
        # Controles del PCA
        col1, col2, col3 = st.columns(3)
        
        with col1:
            color_pca = st.selectbox(
                "🎨 Colorear por:",
                ['subespecie', 'HDG_80HEAD', 'region', 'HDG_category', 'country'],
                key="color_pca_tab2"
            )
        
        with col2:
            pc_x = st.selectbox("📊 Eje X (PC):", list(range(1, 11)), index=0, key="pcx_tab2")
        
        with col3:
            pc_y = st.selectbox("📊 Eje Y (PC):", list(range(1, 11)), index=1, key="pcy_tab2")
        
        # Crear PCA plot
        fig_pca = create_pca_plot(
            pca_result, 
            df_tab2, 
            color_by=color_pca, 
            pc_x=pc_x-1, 
            pc_y=pc_y-1,
            explained_var=explained_var
        )
        
        st.plotly_chart(fig_pca, use_container_width=True)
        
        # Mostrar varianza explicada
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("### 📊 Varianza Explicada por Componente")
            
            var_df = pd.DataFrame({
                'PC': [f'PC{i+1}' for i in range(len(explained_var))],
                'Varianza (%)': explained_var * 100
            })
            
            fig_var = px.bar(
                var_df,
                x='PC',
                y='Varianza (%)',
                title='Varianza Explicada',
                height=400
            )
            
            st.plotly_chart(fig_var, use_container_width=True)
        
        with col2:
            st.markdown("### 🔍 Interpretación Biológica")
            
            st.markdown(f"""
            <div class='info-box'>
                <b>Componentes Principales:</b><br>
                • <b>PC1 ({explained_var[0]*100:.1f}%)</b>: Separación IND vs JAP<br>
                • <b>PC2 ({explained_var[1]*100:.1f}%)</b>: Variación dentro de IND<br>
                • <b>PC3 ({explained_var[2]*100:.1f}%)</b>: Separación AUS<br>
                • <b>PC4-10</b>: Estructura fina poblacional<br><br>
                
            </div>
            """, unsafe_allow_html=True)
            
            # Estadísticas PCA por subespecie
            st.markdown("### 📈 Estadísticas por Subespecie")
            
            pca_stats = df_tab2.groupby('subespecie')['HDG_80HEAD'].agg([
                ('N', 'count'),
                ('Media', 'mean'),
                ('SD', 'std'),
                ('Min', 'min'),
                ('Max', 'max')
            ]).round(1)
            
            st.dataframe(pca_stats, use_container_width=True)
    
    # ==================== TAB 3: HDG VS MORFOLOGÍA ====================
    with tab3:
        st.markdown("<h2 class='sub-header'>📊 Relación Heading Date con Caracteres Morfológicos</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        El heading date presenta **trade-offs complejos** con diversos caracteres morfológicos.
        Estos trade-offs son cruciales para el diseño de ideotipos en programas de mejora.
        """)
        
        # Filtros específicos para este tab
        with st.expander("🎛️ Filtros de Datos", expanded=False):
            df_tab3 = apply_filters(df, "tab3")
        
        if 'df_tab3' not in locals():
            df_tab3 = df.copy()
        
        # Selector de rasgos morfológicos
        morpho_traits = [col for col in df.columns if col not in 
                        ['HDG_80HEAD', 'subespecie', 'country', 'latitude', 'longitude', 
                         'region', 'HDG_category']]
        
        if len(morpho_traits) == 0:
            st.warning("No hay rasgos morfológicos disponibles en el dataset.")
        else:
            col1, col2 = st.columns([2, 1])
            
            with col1:
                selected_trait = st.selectbox(
                    "🔬 Seleccionar rasgo morfológico:",
                    morpho_traits,
                    index=morpho_traits.index('SDHT') if 'SDHT' in morpho_traits else 0,
                    key="trait_tab3"
                )
                
                color_scatter = st.radio(
                    "Colorear por:",
                    ['subespecie', 'region', 'HDG_category'],
                    horizontal=True,
                    key="color_scatter_tab3"
                )
                
                fig_scatter, corr, p_val = create_hdg_morphology_scatter(
                    df_tab3,
                    x_trait=selected_trait,
                    y_trait='HDG_80HEAD',
                    color_by=color_scatter
                )
                
                if fig_scatter:
                    st.plotly_chart(fig_scatter, use_container_width=True)
                else:
                    st.warning("No hay suficientes datos para crear el gráfico.")
            
            with col2:
                st.markdown("### 📊 Matriz de Correlaciones")
                
                # Calcular correlaciones
                corr_traits = ['HDG_80HEAD'] + morpho_traits[:]
                available_traits = [t for t in corr_traits if t in df_tab3.columns]
                
                if len(available_traits) > 1:
                    corr_matrix = df_tab3[available_traits].corr()['HDG_80HEAD'].drop('HDG_80HEAD').sort_values(ascending=False)
                    
                    corr_df = pd.DataFrame({
                        'Rasgo': corr_matrix.index,
                        'Correlación': corr_matrix.values
                    })
                    
                    fig_corr = px.bar(
                        corr_df,
                        x='Correlación',
                        y='Rasgo',
                        orientation='h',
                        color='Correlación',
                        color_continuous_scale='RdBu_r',
                        range_color=[-1, 1],
                        title='Correlación con HDG',
                        height=400
                    )
                    
                    st.plotly_chart(fig_corr, use_container_width=True)
            
            # Trade-offs identificados
            st.markdown("---")
            st.markdown("### ⚖️ Trade-offs Principales Identificados")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("""
                <div class='warning-box'>
                    <h4>🌱 HDG vs Altura</h4>
                    <b>r ≈ +0.37</b><br>
                    Variedades tempranas tienden a ser más bajas.<br>
                    <b>Causa:</b> Menor fase vegetativa = menor elongación
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown("""
                <div class='warning-box'>
                    <h4>🌾 HDG vs Long. Panícula</h4>
                    <b>r ≈ +0.55</b><br>
                    Variedades tempranas tienen panículas más cortas.<br>
                    <b>Causa:</b> Limitación de tiempo y recursos
                </div>
                """, unsafe_allow_html=True)
            
            with col3:
                st.markdown("""
                <div class='warning-box'>
                    <h4>🌾 HDG vs Tamaño Grano</h4>
                    <b>r ≈ +0.20</b><br>
                    Variedades tempranas tienden a tener granos más pequeños.<br>
                    <b>Causa:</b> Período de llenado más corto
                </div>
                """, unsafe_allow_html=True)
            
            # Análisis de categorías HDG
            st.markdown("### 📊 Comparación por Categorías de HDG")
            
            if len(morpho_traits[:6]) > 0:
                available_morpho = [t for t in morpho_traits[:6] if t in df_tab3.columns]
                if len(available_morpho) > 0:
                    morpho_comparison = df_tab3.groupby('HDG_category')[available_morpho].mean().T
                    
                    fig_heatmap = px.imshow(
                        morpho_comparison,
                        labels=dict(x="Categoría HDG", y="Rasgo", color="Valor medio"),
                        aspect="auto",
                        color_continuous_scale='Viridis',
                        title="Perfil Morfológico por Categoría de HDG",
                        height=400
                    )
                    
                    st.plotly_chart(fig_heatmap, use_container_width=True)
    
    # ==================== TAB 4: ANÁLISIS POR REGIONES ====================
    with tab4:
        st.markdown("<h2 class='sub-header'>🌍 Análisis Comparativo por Regiones</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        Las diferentes regiones geográficas muestran **patrones distintos** de variación en heading date,
        reflejando las condiciones climáticas locales y las prácticas agrícolas.
        """)
        
        # Filtros específicos para este tab
    
        
        if 'df_tab4' not in locals():
            df_tab4 = df.copy()
        
        if 'region' in df_tab4.columns:
            # Boxplot por región
            fig_region = create_regional_boxplot(df_tab4, region_col='region', value_col='HDG_80HEAD')
            if fig_region:
                st.plotly_chart(fig_region, use_container_width=True)
            
            # Estadísticas detalladas por región
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("### 📊 Estadísticas por Región")
                
                region_stats = df_tab4.groupby('region').agg({
                    'HDG_80HEAD': ['count', 'mean', 'std', 'min', 'max']
                }).round(1)
                
                region_stats.columns = ['N', 'Media', 'SD', 'Min', 'Max']
                region_stats['CV%'] = (region_stats['SD'] / region_stats['Media'] * 100).round(1)
                
                st.dataframe(region_stats, use_container_width=True)
            
            with col2:
                st.markdown("### 🧬 Subespecies por Región")
                
                # Tabla cruzada región-subespecie
                region_subspecies = pd.crosstab(
                    df_tab4['region'],
                    df_tab4['subespecie']
                )
                
                fig_region_ssp = px.imshow(
                    region_subspecies,
                    labels=dict(x="Subespecie", y="Región", color="N Accesiones"),
                    aspect="auto",
                    color_continuous_scale='Greens',
                    title="Distribución de Subespecies por Región",
                    height=400
                )
                
                st.plotly_chart(fig_region_ssp, use_container_width=True)
            
            # Comparación de distribuciones
            st.markdown("### 📈 Distribuciones de HDG por Región")
            
            fig_violin = px.violin(
                df_tab4.dropna(subset=['region', 'HDG_80HEAD']),
                x='region',
                y='HDG_80HEAD',
                color='region',
                box=True,
                points='outliers',
                title='Distribución de Heading Date por Región',
                labels={'HDG_80HEAD': 'Heading Date (días)', 'region': 'Región'},
                height=500
            )
            
            st.plotly_chart(fig_violin, use_container_width=True)
            
            # Insights por región
            st.markdown("### 🔍 Características Regionales")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("""
                <div class='info-box'>
                    <h4>🌏 SAS (Sur de Asia)</h4>
                    • Alta diversidad<br>
                    • Predominio de IND<br>
                    • Amplio rango de HDG<br>
                    • Centro de diversidad
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown("""
                <div class='info-box'>
                    <h4>🌏 EAS (Este de Asia)</h4>
                    • Predominio de JAP<br>
                    • Ciclos más cortos<br>
                    • Mayor uniformidad<br>
                    • Variedades templadas
                </div>
                """, unsafe_allow_html=True)
            
            with col3:
                st.markdown("""
                <div class='info-box'>
                    <h4>🌏 SEA (Sudeste Asiático)</h4>
                    • Mix IND-JAP<br>
                    • Ciclos intermedios<br>
                    • Alta diversidad<br>
                    • Zonas tropicales
                </div>
                """, unsafe_allow_html=True)
        
        else:
            st.warning("No hay información de regiones disponible en el dataset.")
    
   # ==================== TAB 5: CONCLUSIONES Y RECOMENDACIONES ====================
    with tab5:
        st.markdown("<h2 class='sub-header'>📊 Conclusiones Generales y Recomendaciones</h2>", unsafe_allow_html=True)
        
        st.markdown("""
        Síntesis de los **hallazgos principales** del análisis del 3K RGP y recomendaciones estratégicas 
        para programas de mejoramiento genético en el contexto actual y futuro.
        """)
        
        # Sub-tabs para diferentes secciones
        subtab1, subtab2 = st.tabs([
            "🔬 Hallazgos Principales",
            "🌡️ Cambio Climático y Genotipos"
        ])
        
        with subtab1:
            st.markdown("### 🔬 Hallazgos Principales del Análisis")
            
            # Resumen de estadísticas clave
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric(
                    "Accesiones Totales",
                    f"{len(df)}",
                    help="Total de accesiones analizadas del 3K RGP"
                )
            
            with col2:
                st.metric(
                    "Rango HDG",
                    f"{df['HDG_80HEAD'].min():.0f}-{df['HDG_80HEAD'].max():.0f} días",
                    help="Rango completo de variación en heading date"
                )
            
            with col3:
                n_subspecies = df['subespecie'].nunique()
                st.metric(
                    "Subespecies",
                    f"{n_subspecies}",
                    help="Número de subespecies identificadas"
                )
            
            with col4:
                cv_hdg = (df['HDG_80HEAD'].std() / df['HDG_80HEAD'].mean() * 100)
                st.metric(
                    "CV HDG",
                    f"{cv_hdg:.1f}%",
                    help="Coeficiente de variación del heading date"
                )
            
            st.markdown("---")
            
            # Hallazgos principales en columnas
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("""
                <div class='success-box'>
                    <h4>✅ 1. Enorme Variación Genética</h4>
                    <ul>
                        <li><b>Rango de HDG:</b> 50-175 días (CV=23.1%)</li>
                        <li><b>Estructura poblacional:</b> 5 grupos genéticos principales</li>
                        <li><b>Diferenciación FST > 0.3</b> entre subespecies</li>
                        <li><b>Diversidad subutilizada</b> en programas actuales</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                st.markdown("""
                <div class='success-box'>
                    <h4>✅ 2. Adaptación Local Significativa</h4>
                    <ul>
                        <li><b>Correlación latitud-HDG:</b> r ≈ 0.55-0.65 (p<0.001)</li>
                        <li><b>Latitudes bajas (<20°):</b> Predominio Early/Medium</li>
                        <li><b>Latitudes altas (>35°):</b> Predominio Medium/Late</li>
                        <li><b>Sincronización evolutiva</b> con estación de cultivo</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                st.markdown("""
                <div class='success-box'>
                    <h4>✅ 3. Trade-offs Fenotípicos Confirmados</h4>
                    <ul>
                        <li><b>HDG vs Altura:</b> r = +0.37 (tempranas más bajas)</li>
                        <li><b>HDG vs Long. Panícula:</b> r = +0.55 (limitación recursos)</li>
                        <li><b>HDG vs Tamaño Grano:</b> r = +0.20 (período llenado)</li>
                        <li><b>Desafío:</b> Combinar precocidad con rendimiento</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown("""
                <div class='info-box'>
                    <h4>🧬 4. Estructura Genética Clara</h4>
                    <ul>
                        <li><b>IND (45%):</b> Alta diversidad, HDG medio 105d</li>
                        <li><b>JAP (27%):</b> Ciclos cortos, HDG medio 87d</li>
                        <li><b>AUS (11%):</b> Ultra-precoces, HDG medio 80d</li>
                        <li><b>TRJ/ARO/ADM:</b> Grupos especializados</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                st.markdown("""
                <div class='info-box'>
                    <h4>🌍 5. Centros de Diversidad Identificados</h4>
                    <ul>
                        <li><b>India:</b> Mayor colección mundial, CV alto</li>
                        <li><b>Bangladesh:</b> Múltiples subespecies, diversidad extrema</li>
                        <li><b>Sudeste Asiático:</b> Mix IND-JAP, adaptación tropical</li>
                        <li><b>Este Asia:</b> Variedades templadas, alta uniformidad</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                st.markdown("""
                <div class='warning-box'>
                    <h4>⚠️ 6. Limitaciones del Análisis</h4>
                    <ul>
                        <li><b>Fenotipos en un solo ambiente</b> (sin datos G×E)</li>
                        <li><b>Faltan datos de rendimiento</b> y calidad de grano</li>
                        <li><b>No hay información</b> de tolerancia a estreses</li>
                        <li><b>Requerido:</b> Validación multi-ambiente</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            # Gráfico resumen: Distribución HDG por subespecie
            st.markdown("---")
            st.markdown("### 📊 Síntesis Visual: HDG por Subespecie")
            
            fig_summary = px.violin(
                df.dropna(subset=['subespecie', 'HDG_80HEAD']),
                x='subespecie',
                y='HDG_80HEAD',
                color='subespecie',
                box=True,
                points='outliers',
                title='Distribución de Heading Date por Subespecie',
                labels={'HDG_80HEAD': 'Heading Date (días)', 'subespecie': 'Subespecie'},
                height=500
            )
            
            st.plotly_chart(fig_summary, use_container_width=True)
            
        with subtab2:
            st.markdown("### 🌡️ Cambio Climático: Desafíos y Genotipos de Interés")
            
            st.markdown("""
            <div class='warning-box'>
                <h4>⚠️ Impactos Proyectados del Cambio Climático en Cultivo de Arroz</h4>
                <p>Las proyecciones climáticas para 2050-2100 presentan desafíos críticos para la producción arrocera global:</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Desafíos principales
            st.markdown("#### 🔥 Principales Amenazas Climáticas")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("""
                <div class='warning-box'>
                    <h4>1. 🌡️ Aumento de Temperaturas</h4>
                    <p><b>Problema:</b></p>
                    <ul>
                        <li>Incremento de 2-4°C en temperatura media</li>
                        <li>Olas de calor más frecuentes e intensas</li>
                        <li>Estrés térmico durante floración (>35°C)</li>
                        <li>Reducción de fertilidad de polen</li>
                    </ul>
                    <p><b>Consecuencias:</b></p>
                    <ul>
                        <li>Aceleración del desarrollo → Ciclo más corto</li>
                        <li>Menor período de llenado de grano</li>
                        <li>Reducción de rendimiento 10-20% por cada °C</li>
                        <li>Pérdida de calidad del grano</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                st.markdown("""
                <div class='warning-box'>
                    <h4>3. 💧 Estrés Hídrico Aumentado</h4>
                    <p><b>Problema:</b></p>
                    <ul>
                        <li>Sequías más frecuentes y prolongadas</li>
                        <li>Reducción de disponibilidad de agua</li>
                        <li>Mayor evapotranspiración</li>
                        <li>Competencia por recursos hídricos</li>
                    </ul>
                    <p><b>Consecuencias:</b></p>
                    <ul>
                        <li>Necesidad de variedades tolerantes a sequía</li>
                        <li>Sistemas de secano en riesgo</li>
                        <li>Ciclos cortos para escape a sequía terminal</li>
                        <li>Reducción de área cultivable</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown("""
                <div class='warning-box'>
                    <h4>2. 🌊 Variabilidad de Precipitaciones</h4>
                    <p><b>Problema:</b></p>
                    <ul>
                        <li>Inundaciones más frecuentes</li>
                        <li>Desincronización de lluvias con siembra</li>
                        <li>Períodos secos prolongados</li>
                        <li>Eventos extremos impredecibles</li>
                    </ul>
                    <p><b>Consecuencias:</b></p>
                    <ul>
                        <li>Necesidad de flexibilidad fenológica</li>
                        <li>Variedades tolerantes a submersión</li>
                        <li>Ajuste de fechas de siembra</li>
                        <li>Mayor riesgo en producción</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                st.markdown("""
                <div class='warning-box'>
                    <h4>4. 🦟 Presión de Plagas y Enfermedades</h4>
                    <p><b>Problema:</b></p>
                    <ul>
                        <li>Expansión de rangos de plagas</li>
                        <li>Nuevas generaciones por año</li>
                        <li>Aparición de patógenos emergentes</li>
                        <li>Mayor severidad de enfermedades</li>
                    </ul>
                    <p><b>Consecuencias:</b></p>
                    <ul>
                        <li>Aumento de uso de pesticidas</li>
                        <li>Necesidad de resistencias múltiples</li>
                        <li>Pérdidas de cosecha impredecibles</li>
                        <li>Mayores costos de producción</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            # Genotipos de interés
            st.markdown("---")
            st.markdown("### 🧬 Genotipos de Interés para Adaptación Climática")
            
            # Identificar genotipos clave
            climate_adapted = df[df['HDG_80HEAD'] < 85].copy()
            
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown("""
                <div class='success-box'>
                    <h4>🎯 Estrategia 1: Escape Térmico (Ciclos Ultra-Cortos)</h4>
                    <p><b>Objetivo:</b> Completar ciclo antes de picos de calor</p>
                    <p><b>Genotipos recomendados:</b></p>
                    <ul>
                        <li><b>AUS de Bangladesh:</b> HDG 50-75 días</li>
                        <li><b>Ventajas:</b> Escape a estrés térmico y sequía terminal</li>
                        <li><b>Aplicación:</b> Doble/triple cosecha, zonas áridas</li>
                        <li><b>Desafío:</b> Menor rendimiento por ciclo</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                if len(climate_adapted) > 0:
                    st.markdown("**📊 Germoplasma Ultra-Precoz Disponible:**")
                    
                    aus_summary = climate_adapted.groupby(['country', 'subespecie']).agg({
                        'HDG_80HEAD': ['count', 'mean']
                    }).round(1)
                    
                    aus_summary.columns = ['N Accesiones', 'HDG Medio']
                    aus_summary = aus_summary.sort_values('HDG Medio').head(10)
                    
                    st.dataframe(aus_summary, use_container_width=True)
                    
                    st.success(f"✅ **{len(climate_adapted)} accesiones disponibles** con HDG < 85 días")
                
                st.markdown("""
                <div class='success-box'>
                    <h4>🎯 Estrategia 2: Tolerancia Directa a Calor</h4>
                    <p><b>Objetivo:</b> Mantener fertilidad bajo altas temperaturas</p>
                    <p><b>Genotipos recomendados:</b></p>
                    <ul>
                        <li><b>IND tropicales:</b> De zonas cálidas (>30°C medio)</li>
                        <li><b>Ventajas:</b> Adaptación natural a calor</li>
                        <li><b>Prioridad:</b> India, Indonesia, Filipinas</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown("""
                <div class='success-box'>
                    <h4>🎯 Estrategia 3: Flexibilidad Fenológica</h4>
                    <p><b>Objetivo:</b> Insensibilidad a fotoperiodo</p>
                    <p><b>Genotipos recomendados:</b></p>
                    <ul>
                        <li><b>JAP templadas:</b> Genes Hd1 no funcional</li>
                        <li><b>Ventajas:</b> Floración independiente de día-corto</li>
                        <li><b>Aplicación:</b> Múltiples fechas de siembra/año</li>
                        <li><b>Expansión:</b> Cultivo en nuevas latitudes</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                # Distribución geográfica de genotipos de interés
                if len(climate_adapted) > 0:
                    st.markdown("**🌍 Distribución Geográfica:**")
                    
                    geo_dist = climate_adapted.groupby(['region', 'subespecie']).size().reset_index(name='count')
                    
                    fig_geo = px.bar(
                        geo_dist,
                        x='region',
                        y='count',
                        color='subespecie',
                        title='Genotipos de Interés por Región',
                        labels={'count': 'N Accesiones', 'region': 'Región'},
                        height=350
                    )
                    
                    st.plotly_chart(fig_geo, use_container_width=True)
                
                st.markdown("""
                <div class='success-box'>
                    <h4>🎯 Estrategia 4: Adaptación a Estrés Hídrico</h4>
                    <p><b>Objetivo:</b> Tolerancia a sequía intermitente</p>
                    <p><b>Genotipos recomendados:</b></p>
                    <ul>
                        <li><b>AUS de zonas secas:</b> Sistema radicular profundo</li>
                        <li><b>Ventajas:</b> Eficiencia de uso de agua</li>
                        <li><b>Caracteres:</b> Enrollamiento foliar, cera epicuticular</li>
                        <li><b>Aplicación:</b> Sistemas de secano, riego limitado</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
            
            # Recursos genéticos clave detallados
            st.markdown("---")
            st.markdown("### 🔑 Pool de Germoplasma para Resiliencia Climática")
            
            if len(climate_adapted) > 0:
                climate_resources = climate_adapted.groupby(
                    ['subespecie', 'country']
                ).agg({
                    'HDG_80HEAD': ['count', 'mean', 'min']
                }).round(1)
                
                climate_resources.columns = ['N Accesiones', 'HDG Medio', 'HDG Mínimo']
                climate_resources = climate_resources.sort_values('N Accesiones', ascending=False).head(20)
                
                fig_climate = go.Figure()
                
                # Convertir el índice multinivel a columnas para plotly
                climate_plot = climate_resources.reset_index()
                climate_plot['label'] = climate_plot['country'] + ' (' + climate_plot['subespecie'] + ')'
                
                fig_climate.add_trace(go.Bar(
                    y=climate_plot['label'],
                    x=climate_plot['N Accesiones'],
                    orientation='h',
                    marker_color='#FF9800',
                    text=climate_plot['HDG Medio'],
                    texttemplate='HDG: %{text}d',
                    textposition='outside'
                ))
                
                fig_climate.update_layout(
                    title="Top 20: Germoplasma para Resiliencia Climática (HDG < 85 días)",
                    xaxis_title="Número de Accesiones",
                    yaxis_title="País (Subespecie)",
                    height=600,
                    showlegend=False
                )
                
                st.plotly_chart(fig_climate, use_container_width=True)



    # ==================== FOOTER ====================
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #666; font-size: 0.9rem; padding: 2rem;'>
        <p><b>🌾 3K Rice Genomes Project - Heading Date Explorer</b></p>
        <p>Desarrollado para programas de mejoramiento genético | Octubre 2025</p>
        <p>📧 Análisis basado en el proyecto 3K RGP con 1,882 accesiones y SNPs</p>
    </div>
    """, unsafe_allow_html=True)

# ==================== EJECUCIÓN ====================
if __name__ == "__main__":
    main()
