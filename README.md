# genomeRice

# 🌾 3K Rice Genomes - Heading Date Explorer

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.28%2B-FF4B4B.svg)](https://streamlit.io/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Data](https://img.shields.io/badge/Data-3K%20RGP-orange.svg)](http://snp-seek.irri.org/)

Aplicación interactiva de análisis exploratorio para visualizar y analizar la variación del **Heading Date (HDG)** y su relación con la estructura genética, geografía y morfología en el arroz asiático cultivado (*Oryza sativa*), utilizando datos del proyecto **3,000 Rice Genomes (3K RGP)**. 
**Dataset : https://www.kaggle.com/datasets/saurabhshahane/rice-genotype**


---

## 📋 Tabla de Contenidos

- [Descripción del Proyecto](#-descripción-del-proyecto)
- [Características Principales](#-características-principales)
- [Instalación](#-instalación)
- [Uso](#-uso)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Dataset](#-dataset)
- [Metodología](#-metodología)
- [Resultados Clave](#-resultados-clave)
- [Visualizaciones](#-visualizaciones)
- [Tecnologías](#-tecnologías)
- [Referencias Científicas](#-referencias-científicas)

---

## 🔬 Descripción del Proyecto

Este proyecto desarrolla una **aplicación web interactiva** construida con Streamlit para explorar y analizar datos del **3,000 Rice Genomes Project (3K RGP)**, que contiene información genómica y fenotípica de **3,010 accesiones** de arroz asiático de todo el mundo.

### Objetivos Principales

1. **Caracterizar la variación del Heading Date (HDG_80HEAD)** como rasgo clave para adaptación climática
2. **Explorar la estructura genética poblacional** mediante análisis de componentes principales (PCA)
3. **Identificar relaciones entre HDG y caracteres morfológicos** (trade-offs fenotípicos)
4. **Analizar patrones de adaptación geográfica** y correlación con latitud
5. **Identificar germoplasma elite** para diferentes objetivos de mejoramiento genético
6. **Proporcionar recomendaciones estratégicas** frente al cambio climático

### ¿Por qué es importante?

El **heading date** (tiempo a floración) es uno de los caracteres más críticos en arroz porque:
- Determina la **adaptación a diferentes zonas agroecológicas**
- Afecta directamente el **rendimiento y la calidad del grano**
- Es clave para **escapar a estreses** (sequía terminal, calor extremo)
- Permite **múltiples cosechas por año** en zonas tropicales
- Es esencial para la **adaptación al cambio climático**

---

## ✨ Características Principales

### 🗺️ **Exploración Geográfica**
- Mapas interactivos con distribución global de accesiones
- Coloración dinámica por HDG, subespecie, región o categoría
- Análisis de correlación latitud-HDG 
- Identificación de centros de diversidad

### 🧬 **Análisis Genético (PCA)**
- PCA interactivo de 12,486 SNPs
- Visualización de estructura poblacional (5 subespecies)
- Selección de componentes principales (PC1-PC10)
- Coloración por múltiples variables (subespecie, HDG, región, país)
- Varianza explicada por cada componente

### 📊 **HDG vs Morfología**
- Gráficos de dispersión interactivos
- Análisis de correlaciones (matriz visual)
- Identificación de trade-offs fenotípicos:
  - HDG vs Altura (r ≈ +0.37)
  - HDG vs Longitud de Panícula (r ≈ +0.55)
  - HDG vs Tamaño de Grano (r ≈ +0.20)
- Comparación por categorías de HDG (Early, Medium, Late)

### 🌍 **Análisis por Regiones**
- Boxplots comparativos por región (SAS, EAS, SEA, AFR, EUR, etc.)
- Distribución de subespecies por región (heatmap)
- Violin plots para visualizar distribuciones completas
- Estadísticas detalladas (media, SD, CV%, rango)

### 🎯 **Conclusiones y Recomendaciones**
- Síntesis de hallazgos principales
- **Análisis exhaustivo del cambio climático**:
  - 4 amenazas principales (temperatura, precipitaciones, sequía, plagas)
  - 4 estrategias de adaptación con genotipos específicos
- Recomendaciones estratégicas por plazos:
  - Corto plazo (1-3 años): MAS, bancos de germoplasma
  - Mediano plazo (3-5 años): Selección genómica, NILs
  - Largo plazo (5-10 años): CRISPR, multi-ómica, breeding digital

---

## 🚀 Instalación

### Requisitos Previos

- Python 3.8 o superior
- pip (gestor de paquetes de Python)
- Git (opcional, para clonar el repositorio)

### Paso 1: Clonar el Repositorio

```bash
git clone https://github.com/tu-usuario/3k-rice-hdg-explorer.git
cd 3k-rice-hdg-explorer
```

### Paso 2: Crear Entorno Virtual (Recomendado)

```bash
# En Windows
python -m venv venv
venv\Scripts\activate

# En macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

### Paso 3: Instalar Dependencias

```bash
pip install -r requirements.txt
```

### Contenido de `requirements.txt`:

```txt
streamlit>=1.28.0
pandas>=2.0.0
numpy>=1.24.0
plotly>=5.17.0
scikit-learn>=1.3.0
scipy>=1.11.0
```

### Paso 4: Preparar los Datos

Coloca tus archivos de datos en la carpeta `reports/`:

```
project/
├── reports/
│   ├── traits_geno_aligned.csv    # Datos fenotípicos (1882 x N)
│   └── geno_aligned.csv            # Datos genotípicos (1882 x 12486)
├── app.py
└── requirements.txt
```

**Nota:** Si no tienes los datos reales, la aplicación generará automáticamente datos de demostración realistas basados en las estadísticas del 3K RGP.

---

## 💻 Uso

### Ejecutar la Aplicación

```bash
streamlit run app.py
```

La aplicación se abrirá automáticamente en tu navegador en `http://localhost:8501`

### Navegación por la Aplicación

1. **🗺️ Exploración Geográfica**: 
   - Expande "🎛️ Filtros de Datos" para aplicar filtros
   - Selecciona el color del mapa (HDG, subespecie, región, categoría)
   - Explora los gráficos de distribución por país
   - Analiza la correlación latitud-HDG

2. **🧬 Análisis Genético (PCA)**:
   - Aplica filtros independientes si es necesario
   - Selecciona la variable de coloración
   - Elige los ejes del PCA (PC1-PC10)
   - Observa la varianza explicada

3. **📊 HDG vs Morfología**:
   - Filtra los datos según tu interés
   - Selecciona el rasgo morfológico a analizar
   - Colorea por subespecie, región o categoría
   - Revisa la matriz de correlaciones

4. **🌍 Análisis por Regiones**:
   - Compara distribuciones entre regiones
   - Analiza la distribución de subespecies
   - Identifica características regionales

5. **🎯 Conclusiones y Recomendaciones**:
   - Explora genotipos de interés para cambio climático
   - Consulta las recomendaciones estratégicas

---

## 📁 Estructura del Proyecto

```
3k-rice-hdg-explorer/
│
├── app.py                          # Aplicación principal de Streamlit
├── requirements.txt                # Dependencias de Python
├── README.md                       # Este archivo
├── LICENSE                         # Licencia del proyecto
│
├── reports/                        # Datos procesados
│   ├── traits_geno_aligned.csv    # Fenotipos alineados
│   └── geno_aligned.csv            # Genotipos alineados
│
├── notebooks/                      # Notebooks de análisis
│   ├── eda_improved.ipynb          # EDA completo
│   └── preprocessing.ipynb         # Preprocesamiento de datos
│
├── docs/                           # Documentación adicional
│   ├── images/                     # Imágenes para README
│   ├── methodology.md              # Metodología detallada
│   └── references.md               # Referencias completas
│
└── data/                           # Datos brutos (no incluidos en repo)
    ├── raw/                        # Datos originales del 3K RGP
    └── processed/                  # Datos intermedios
```

---

## 📊 Dataset

### 3,000 Rice Genomes Project (3K RGP)

El **3K RGP** es el proyecto de secuenciación de arroz más grande hasta la fecha, coordinado por el International Rice Research Institute (IRRI) y publicado en **Nature (2018)**.

#### Características del Dataset:

- **Accesiones**: 1,882 variedades de arroz asiático (*Oryza sativa*)
- **Origen**: 89 países de Asia, África, América y Europa
- **SNPs**: 12,486 marcadores moleculares (filtrados de 29 millones)
- **Cobertura genómica**: ~40x por muestra
- **Subespecies principales**:
  - **IND** (Indica): 45% - zonas tropicales
  - **JAP** (Japonica): 27% - zonas templadas
  - **AUS** (Aus): 11% - Bangladesh, India
  - **ARO** (Aromatic): 4% - Basmati, Jasmine
  - **TRJ** (Tropical Japonica): 9% - África, América Latina
  - **ADM** (Admixed): 4% - híbridos

#### Variables Fenotípicas Principales:

| Variable | Descripción | Unidad | Rango |
|----------|-------------|--------|-------|
| `HDG_80HEAD` | Días a 80% de floración | días | 50-175 |
| `SDHT` | Altura de plántula | cm | 20-60 |
| `PLT_POST` | Longitud de panícula | cm | 15-35 |
| `CULT_REPRO` | Duración fase reproductiva | días | 60-150 |
| `GRLT` | Longitud de grano | mm | 4-12 |
| `GRWD` | Ancho de grano | mm | 1.5-4 |
| `GRWT100` | Peso de 100 granos | g | 1-4 |
| `LLT` | Longitud de hoja bandera | cm | 25-70 |
| `LWD` | Ancho de hoja bandera | cm | 0.8-2.5 |

#### Fuentes de Datos:


- **Portal oficial**: [SNP-Seek Database](http://snp-seek.irri.org/)
- **Publicación principal**: Wang et al. (2018) *Nature* 557:43-49
- **Datos genómicos**: [NCBI BioProject PRJNA301661](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301661)

---

## 🧪 Metodología

### 1. Preprocesamiento de Datos

```python
# Alineación de fenotipos y genotipos
- Eliminación de valores faltantes
- Normalización de nombres de países
- Cálculo de variables derivadas (GrainSize, HDG_category)
- Geocodificación de países
```

### 2. Análisis Exploratorio

- **Estadísticas descriptivas**: Media, mediana, SD, CV, rango
- **Distribuciones**: Histogramas, boxplots, violin plots
- **Correlaciones**: Pearson, Spearman (HDG vs rasgos morfológicos)
- **Análisis geográfico**: Correlación latitud-HDG

### 3. Análisis de Estructura Genética

```python
# PCA en 12,486 SNPs
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# 1. Normalización
scaler = StandardScaler()
geno_scaled = scaler.fit_transform(geno)

# 2. PCA
pca = PCA(n_components=10)
pca_result = pca.fit_transform(geno_scaled)

# 3. Varianza explicada
explained_var = pca.explained_variance_ratio_
```

### 4. Visualizaciones Interactivas

- **Mapas geográficos**: Plotly Mapbox (OpenStreetMap)
- **Gráficos de dispersión**: Plotly Express con trendlines (OLS)
- **Heatmaps**: Correlaciones, distribuciones cruzadas
- **Boxplots/Violin plots**: Comparaciones entre grupos



---

## 🔑 Resultados Clave

### Variación Genética

- **Rango de HDG**: 50-175 días (CV = 23.1%)
- **Estructura poblacional**: 5 grupos con FST > 0.3
- **Diversidad**: Subutilizada en programas actuales


### Centros de Diversidad

1. **India**: Mayor colección mundial (>500 accesiones), CV alto
2. **Bangladesh**: Múltiples subespecies, diversidad extrema (AUS)
3. **Indonesia**: Diversidad tropical única (IND + JAP)
4. **China**: Amplio rango geográfico, variación latitudinal


---

## 🛠️ Tecnologías

### Lenguajes y Frameworks

- **Python 3.8+**: Lenguaje principal
- **Streamlit 1.28+**: Framework web interactivo
- **Plotly 5.17+**: Visualizaciones interactivas

### Librerías de Análisis

- **Pandas**: Manipulación de datos
- **NumPy**: Operaciones numéricas
- **Scikit-learn**: Machine Learning (PCA, StandardScaler)
- **SciPy**: Estadística (correlaciones, tests)

### Visualización

- **Plotly Express**: Gráficos rápidos
- **Plotly Graph Objects**: Gráficos personalizados
- **Mapbox**: Mapas geográficos interactivos

---

## 🤝 Contribuciones

¡Las contribuciones son bienvenidas! Por favor, sigue estos pasos:

1. **Fork** el repositorio
2. Crea una **rama** para tu feature (`git checkout -b feature/NuevaCaracteristica`)
3. **Commit** tus cambios (`git commit -m 'Agrega nueva característica'`)
4. **Push** a la rama (`git push origin feature/NuevaCaracteristica`)
5. Abre un **Pull Request**

### Áreas de Mejora

- [ ] Integrar datos de rendimiento y calidad de grano
- [ ] Añadir análisis de interacción genotipo × ambiente (G×E)
- [ ] Implementar modelos de selección genómica
- [ ] Agregar predicciones de fenotipos con ML
- [ ] Incluir análisis de haplotipos en genes candidatos
- [ ] Desarrollar módulo de diseño de cruces
- [ ] Añadir exportación de resultados (PDF, Excel)
- [ ] Internacionalización (i18n) - soporte multiidioma

---

---

```

---


---

## 🙏 Agradecimientos

- **International Rice Research Institute (IRRI)** por liderar el proyecto 3K RGP
- **The 3,000 Rice Genomes Consortium** por hacer los datos públicos
- **Comunidad de mejoramiento genético de arroz** por el valioso feedback
- **Streamlit** por la excelente plataforma de desarrollo
- **Plotly** por las herramientas de visualización interactiva

---



---

<div align="center">

### 🌾 Hecho con ❤️ para la comunidad científica del arroz

**El futuro de la producción arrocera depende de nuestra capacidad de actuar ahora**

</div>
