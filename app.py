import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem, Draw
from rdkit.Chem.Draw import MolToImage
import pandas as pd
from io import BytesIO
import py3Dmol
import streamlit.components.v1 as components

# ─────────────────────────────────────────────────────────────
# 🎨 INSTAGRAM-STYLE CSS & CONFIG
# ─────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="MolDesign Studio",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .main {
        background-color: #fafafa;
    }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    h1 {
        color: #262626;
        font-weight: 700;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        text-align: center;
        margin-bottom: 1rem;
    }
    h2, h3 {
        color: #262626;
        font-weight: 600;
    }
    .css-1r6slb0 {
        background-color: white;
        border-radius: 15px;
        padding: 20px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.05);
        border: 1px solid #efefef;
    }
    div.stButton > button {
        background: linear-gradient(45deg, #f09433 0%, #e6683c 25%, #dc2743 50%, #cc2366 75%, #bc1888 100%); 
        color: white;
        border: none;
        border-radius: 8px;
        padding: 10px 24px;
        font-weight: 600;
        transition: all 0.3s ease;
        width: 100%;
    }
    div.stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(220, 39, 67, 0.4);
    }
    input[type="text"] {
    color: #000000 !important;
    border-radius: 8px;
    border: 1px solid #dbdbdb;
    padding: 10px;
    background-color: #fafafa;
}
    input[type="text"]:focus {
        border-color: #a8a8a8;
        outline: none;
    }
    [data-testid="stMetricValue"] {
        font-size: 1.5rem;
        color: #FFFFFF !important;
    }
    [data-testid="stMetricLabel"] {
        font-size: 0.9rem;
        color: #FFFFFF !important;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    .dataframe {
        border: none !important;
        border-radius: 10px;
        overflow: hidden;
    }
    section[data-testid="stSidebar"] {
        background-color: #ffffff;
        border-right: 1px solid #efefef;
        /* Красный цвет для заголовков метрик (MW, LogP, TPSA, Lipinski) */
[data-testid="stMetricLabel"] {
    color: #e63946 !important;
    font-weight: 600 !important;
}

/* Красный цвет для значений метрик (числа, PASS и т.д.) */
[data-testid="stMetricValue"] {
    color: #e63946 !important;
    font-weight: 700 !important;
}

/* Дополнительно — для дельты (если она есть) */
[data-testid="stMetricDelta"] {
    color: #e63946 !important;
}
    }
</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────
# ⚙️ ФУНКЦИИ
# ─────────────────────────────────────────────────────────────

@st.cache_data
def calculate_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid Structure"}
        
        Chem.SanitizeMol(mol)
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot = Descriptors.NumRotatableBonds(mol)
        
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        
        return {
            "mol": mol,
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "TPSA": round(tpsa, 2),
            "HBD": hbd,
            "HBA": hba,
            "RotBonds": rot,
            "Violations": violations,
            "InChIKey": Chem.MolToInchiKey(mol)
        }
    except Exception as e:
        return {"error": str(e)}

def generate_3d_html(mol, width=600, height=400):
    mol_3d = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol_3d)
        mb = Chem.MolToMolBlock(mol_3d)
    except:
        return "<div style='display:flex;justify-content:center;align-items:center;height:100%;color:#cc2366'>⚠️ 3D Error</div>"
    
    view = py3Dmol.view(width=width, height=height)
    view.addModel(mb, 'mol')
    view.setStyle({'stick': {'radius': 0.15}})
    view.zoomTo()
    return view._make_html()

# ─────────────────────────────────────────────────────────────
# 📱 ИНТЕРФЕЙС
# ─────────────────────────────────────────────────────────────

st.title("🧪 Molecular Design Studio")
st.markdown("<p style='text-align:center; color:#8e8e8e; margin-top:-10px;'>Draw • Analyze • Export</p>", unsafe_allow_html=True)

with st.sidebar:
    st.markdown("<h2 style='color:#000000;'>✨ New Search</h2>", unsafe_allow_html=True)
    
    examples = {
  
        "Аспирин": "CC(=O)Oc1ccccc1C(=O)O",
        "Ибупрофен": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Парацетамол": "CC(=O)Nc1ccc(O)cc1",
        "Диклофенак": "OC(=O)CC1=CC=CC=C1NC2=C(Cl)C=C(Cl)C=C2",
        "Напроксен": "COc1ccc2cc(C(C)C(=O)O)on2c1",
        "Кетопрофен": "CC(C1=CC=CC=C1)C(=O)CC2=CC=CC=C2C(=O)O",
        "Индометацин": "CC1=C(C2=C(C1=O)NC3=C(C2)C=CC=C3Cl)C(=O)OC",
        "Пироксикам": "CN1C(=O)C2=C(C1=O)C=CC=C2S(=O)(=O)N",
        "Целекоксиб": "CS(=O)(=O)C1=CC=C(C=C1)C2=NN(C(=O)C2)C3=CC=C(C=C3)C(F)(F)F",
        "Мелоксикам": "CN1C(=O)C2=C(C1=O)C=CC=C2S(=O)(=O)NC3=C(C(N(C3)=O)=C)C",
        "Пенициллин G": "CC1(C2C(N(C2=O)C1C(=O)O)C(=O)NCc3ccccc3)S",
        "Амоксициллин": "CC(C)(S[C@@H]1[C@H](NC(=O)[C@H](N)C2=CC=C(O)C=C2)C(=O)N1C)C(=O)O",
        "Ампициллин": "CC(C)(S[C@@H]1[C@H](NC(=O)[C@H](N)C2=CC=CC=C2)C(=O)N1C)C(=O)O",
        "Тетрациклин": "CN(C)C1=C(O)C2=C(C(=O)C3=C(O)C4=C(C(O)=C(C(N)=O)C(O)=C4C)C(O)=C3C(=O)N2)C(O)=C1",
        "Эритромицин": "CCC(C)C(OC(=O)CC(O)C(C)C(OC(=O)C(C)C(O)C(C)C(=O)O)CC)OC1OC(C)C(O)C(OC)C1OC2C(C)(O)C(C)C(OC)C2",
        "Ципрофлоксацин": "OC(=O)C1=CC2=C(N3CCNCC3)C(=O)C(C(=O)O)=CN2C=C1",
        "Левофлоксацин": "CN1C2=C(C(=O)C(C(=O)O)=CN1C3CCN(C)CC3)C=C(C=C2)F",
        "Азитромицин": "CN(C)C1(C)OC2CC(C)OC3(C)OC(C)C(O)C(C)C(=O)O)C(O)C(C)C(=O)O)C(O)C(C)C(=O)O)C2",
        "Ванкомицин": "C[C@H]1NC(=O)[C@@H](C(C)C)NC(=O)[C@@H](Cc2ccccc2)NC(=O)[C@H](Cc3ccc(O)cc3)NC(=O)[C@H](Cc4ccc(O)c(O)c4)NC(=O)[C@H](CO)NC(=O)[C@H](Cc5ccc(O)c(O)c5)NC1=O",
        "Атенолол": "CC(C)NC[C@H](O)COc1ccc(cc1)CC(=O)N",
        "Метопролол": "COc1ccc(cc1)C[C@@H](O)CNC(C)C",
        "Пропранолол": "CC(C)NC[C@H](O)COc1cccc2ccccc12",
        "Лизиноприл": "NCCCC[C@H](N)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
        "Эналаприл": "CCOC(=O)[C@H](C)NC(=O)[C@@H](N)CC1=CC=CC=C1",
        "Каптоприл": "SC[C@@H](C)C(=O)N1CCC[C@H]1C(=O)O",
        "Варфарин": "CC(=O)CC(C1=CC=CC=C1)C2=C(O)C3=CC=CC=C3C2=O",
        "Дигоксин": "C[C@H]1OC(O)C(O)C(O)C1OC2CCC3(C)C4(C)CCC5C(C(C5C4CCC3(C)C2)OC6C(C(C(C(O6)C)O)O)O)OC7C(C(C(C(O7)C)O)O)O",
        "Нифедипин": "COC(=O)C1=C(C)NC(=C(C1C2=CC=CC=C2[N+](=O)[O-])C)C3=CC=CC=C3",
        "Верапамил": "COc1ccc(cc1)C(C#N)C(C)N(C)CCC(C)c2cc(OC)c(OC)c(OC)c2",
        "Флуоксетин": "CNCCC(Oc1ccc(cc1)C(F)(F)F)c2ccccc2",
        "Сертралин": "CN[C@H]1CC[C@H]2c3ccc(Cl)cc3CC[C@]21C",
        "Пароксетин": "FC(F)(F)Oc1ccc(cc1)C2CNCCC2c3ccc(O)c(O)c3",
        "Диазепам": "CN1C(=O)CN=C(C2=CC=CC=C2)C3=C1C=C(Cl)C=C3",
        "Лоразепам": "OC1=CC=CC=C1N2C(=O)CN=C(C3=CC=CC=C3)C4=C2C=C(Cl)C=C4",
        "Алпразолам": "CN1C(=O)CN=C(C2=CC=CC=C2)C3=C1C=CC=C3N4N=CC=C4",
        "Клоназепам": "[O-][N+](=O)C1=CC=C(C=C1)N2C(=O)CN=C(C3=CC=CC=C3)C4=C2C=C(Cl)C=C4",
        "Галоперидол": "OC(c1ccc(cc1)C(=O)CCC2CCN(CC2)c3ccc(Cl)cc3Cl)(c4ccc(F)cc4)",
        "Хлорпромазин": "CN(C)CCC1=CC2=C(C=C1)N(C)C3=C2C=CC=C3Cl",
        "Литий карбонат": "[Li+].[Li+].[O-]C(=O)[O-]",
        "Инсулин (фрагмент)": "CC(C)S[C@@H]1C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CS)C(=O)O",
        "Тестостерон": "C[C@@]12CC[C@H](O)C[C@@H]1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CC[C@H](C4)=O)C",
        "Эстрадиол": "C[C@]12CC[C@H](O)C[C@@H]1CCC3=C2CCC4=C3CC[C@@H](O)C4",
        "Прогестерон": "CC(=O)[C@@H]1CC[C@@H]2[C@@H]1CC[C@@H]3[C@@H]2CCC4=CC(=O)CC[C@]34C",
        "Кортизол": "CC12CCC(=O)C=C1CC3C(C2CCC3(CO)O)O",
        "Альдостерон": "CC12CCC(=O)C=C1CC3C(C2CCC3(C=O)O)O",
        "Тироксин": "OC(=O)C(N)Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1",
        "Окситоцин": "N[C@@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CO)C(=O)O1",
        "Вазопрессин": "N[C@@H]1C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CO)C(=O)O1",
        "Кофеин": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Никотин": "CN1CCC[C@H]1c2cccnc2",
        "Морфин": "CN1CCC23C4C1CC(C5=C2C(=C(C=C5)O)O3)O4",
        "Кодеин": "COc1cc2C3CCN(C)CC3Cc4c(O)c(OC)c5c(c12)OCO5",
        "Героин": "CC(=O)Oc1cc2C3CCN(C)CC3Cc4c(O)c(OC)c5c(c12)OCO5",
        "Хинин": "COC1=CC2=C(C=C1)C3CC4(C2)C(C3)N5CCC(C5)C=C6C7=CN=C6C=C7",
        "Хинидин": "COC1=CC2=C(C=C1)C3CC4(C2)C(C3)N5CCC(C5)C=C6C7=CN=C6C=C7",
        "Кокаин": "CN1C2CCC1CC(C2)C(=O)OC3=CC=CC=C3",
        "Атропин": "CN1C2CCC1CC(C2)OC(=O)C(CO)C3=CC=CC=C3",
        "Скополамин": "CN1C2CCC1CC(C2)OC(=O)C(CO)C3=CC=CC=C3",
        "Эфедрин": "C[C@H](N[C@H](C)CO)c1ccccc1",
        "Псевдоэфедрин": "C[C@@H](N[C@H](C)CO)c1ccccc1",
        "Мескалин": "COc1cc(CCN)cc(OC)c1OC",
        "Псилоцибин": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "ЛСД": "CN1CCC2(CC1)c3c[nH]c4cccc3C2N5C(=O)C6=CC=CC=C6C5",
        "Стрихнин": "CN1CCC2(CC1)c3c[nH]c4cccc3C2N5C(=O)C6=CC=CC=C6C5",
        "Бруцин": "COC1=CC2=C(C=C1)C3CC4(C2)C(C3)N5CCC(C5)C=C6C7=CN=C6C=C7",
        "Лимонен": "CC1=CCC(CC1)C(=C)C",
        "Пинен": "CC1=CCC2CC1C2(C)C",
        "Камфора": "CC1(C)C2CCC1(C)C(=O)C2",
        "Ментол": "CC(C)C1CCC(C)C(O)C1",
        "Цитраль": "CC(C)=CCCC(C)=CC=O",
        "Линалоол": "CC(C)=CCCC(C)(O)C=C",
        "Гераниол": "CC(C)=CCCC(C)=CCO",
        "Нерол": "CC(C)=CCCC(C)=CCO",
        "Цитронеллол": "CC(C)CCCC(C)CCO",
        "Эвкалиптол": "CC1CCC2(C)OC1CC2(C)C",
        "Туйон": "CC1CC2CC1C(=O)C2(C)C",
        "Карвакрол": "CC(C)C1=CC(C)=C(O)C=C1",
        "Тимол": "CC(C)C1=CC(C)=C(O)C=C1",
        "Эвгенол": "COc1ccc(CC=C)cc1O",
        "Анетол": "COc1ccc(CC=C)cc1",
        "Кверцетин": "OC1=C(C2=CC=C(O)C=C2OC3=C(O)C=C(O)C=C3O)C(=O)C4=C(O)C=C(O)C=C4O1",
        "Рутин": "OC1=C(C2=CC=C(O)C=C2OC3=C(O)C=C(O)C=C3O)C(=O)C4=C(O)C=C(O)C=C4O1",
        "Катехин": "OC1C(O)C2=CC(O)=C(O)C=C2OC1C3=CC=C(O)C(O)=C3",
        "Эпикатехин": "OC1C(O)C2=CC(O)=C(O)C=C2OC1C3=CC=C(O)C(O)=C3",
        "Апигенин": "OC1=CC(O)=C(C2=CC(=O)C3=C(O)C=C(O)C=C3O2)C=C1",
        "Лютеолин": "OC1=C(C2=CC=C(O)C=C2OC3=C(O)C=C(O)C=C3O)C(=O)C4=C(O)C=C(O)C=C4O1",
        "Нарингенин": "OC1=CC(O)=C(C2=CC(=O)C3=C(O)C=C(O)C=C3O2)C=C1",
        "Гесперетин": "COC1=C(O)C=C(C2=CC(=O)C3=C(O)C=C(O)C=C3O2)C=C1",
        "Изорамнетин": "COC1=C(O)C=C(C2=C(O)C(=O)C3=C(O)C=C(O)C=C3O2)C=C1",
        "Мирицетин": "OC1=C(C2=CC=C(O)C=C2OC3=C(O)C=C(O)C=C3O)C(=O)C4=C(O)C=C(O)C=C4O1",
        "Глюкоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Фруктоза": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)C1=O",
        "Галактоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Манноза": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
        "Рибоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
        "Дезоксирибоза": "OC[C@H]1O[C@H](O)[C@@H](O)C1",
        "Ксилоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
        "Арабиноза": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H]1O",
        "Ликсоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
        "Аллоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Сахароза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2(CO)O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]2O",
        "Лактоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2[C@@H](O)[C@@H](O)[C@H](O)O[C@H]2CO",
        "Мальтоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO",
        "Целлобиоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO",
        "Трегалоза": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO",
        "Аланин": "C[C@H](N)C(=O)O",
        "Аргинин": "NC(N)=NCCC[C@H](N)C(=O)O",
        "Аспарагин": "NC(=O)C[C@H](N)C(=O)O",
        "Аспарагиновая кислота": "OC(=O)C[C@H](N)C(=O)O",
        "Валин": "CC(C)[C@H](N)C(=O)O",
        "Гистидин": "OC(=O)[C@@H](N)Cc1cnc[nH]1",
        "Глицин": "NCC(O)=O",
        "Глутамин": "NC(=O)CC[C@H](N)C(=O)O",
        "Глутаминовая кислота": "OC(=O)CC[C@H](N)C(=O)O",
        "Изолейцин": "CC[C@H](C)[C@H](N)C(=O)O",
        "Лейцин": "CC(C)C[C@H](N)C(=O)O",
        "Лизин": "NCCCC[C@H](N)C(=O)O",
        "Метионин": "CSCC[C@H](N)C(=O)O",
        "Пролин": "OC(=O)[C@@H]1CCCN1",
        "Серин": "OC[C@H](N)C(=O)O",
        "Тирозин": "OC(=O)[C@@H](N)Cc1ccc(O)cc1",
        "Треонин": "C[C@H](O)[C@H](N)C(=O)O",
        "Триптофан": "OC(=O)[C@@H](N)Cc1c[nH]c2ccccc12",
        "Фенилаланин": "OC(=O)[C@@H](N)Cc1ccccc1",
        "Цистеин": "SC[C@H](N)C(=O)O",
        "ГАМК (GABA)": "NCCC(O)=O",
        "Таурин": "NCCS(O)(=O)=O",
        "Орнитин": "NCCCC(N)C(=O)O",
        "Цитруллин": "NC(=O)NCCC(N)C(=O)O",
        "Бета-аланин": "NCCC(O)=O",
        "Креатин": "CN(C)C(=N)NCC(O)=O",
        "Карнитин": "C[C@H](CC(=O)O)C[N+](C)(C)C",
        "Гидроксипролин": "OC1CCNC1C(=O)O",
        "Селеноцистеин": "[Se]C[C@H](N)C(=O)O",
        "Пирролизин": "CC1CCCN1C(=O)C(N)CC(=O)O",
        "Пальмитиновая кислота": "CCCCCCCCCCCCCCCC(=O)O",
        "Стеариновая кислота": "CCCCCCCCCCCCCCCCCC(=O)O",
        "Олеиновая кислота": "CCCCCCCC/C=C\CCCCCCCC(=O)O",
        "Линолевая кислота": "CCCCC/C=C\C/C=C\CCCCCCCC(=O)O",
        "Линоленовая кислота": "CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O",
        "Арахидоновая кислота": "CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC(=O)O",
        "Эйкозапентаеновая кислота (EPA)": "CC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CCCC(=O)O",
        "Докозагексаеновая кислота (DHA)": "CCC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CCC(=O)O",
        "Миристиновая кислота": "CCCCCCCCCCCCCC(=O)O",
        "Лауриновая кислота": "CCCCCCCCCCCC(=O)O",
        "Каприловая кислота": "CCCCCCCC(=O)O",
        "Каприновая кислота": "CCCCCCCCCC(=O)O",
        "Масляная кислота": "CCCC(=O)O",
        "Пропионовая кислота": "CCC(=O)O",
        "Уксусная кислота": "CC(=O)O",
        "Муравьиная кислота": "C(=O)O",
        "Триолеин": "CCCCCCCC/C=C\CCCCCCCC(=O)OCC(OC(=O)CCCCCCC/C=C\CCCCCCCC)COC(=O)CCCCCCC/C=C\CCCCCCCC",
        "Трипальмитин": "CCCCCCCCCCCCCCCC(=O)OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC",
        "Тристеарин": "CCCCCCCCCCCCCCCCCC(=O)OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCC",
        "Холестерин": "CC(C)CCCC(C)C1CCC2C1=CCC3C4CCC(C)(C)CC4CCC23",
        "Холестанол": "CC(C)CCCC(C)C1CCC2C1CCC3C4CCC(C)(C)CC4CCC23",
        "Эргостерин": "CC(C)C=CC(C)C1CCC2C1=CCC3C4CCC(C)(C)CC4CCC23",
        "Ланолин": "CC(C)CCCC(C)C1CCC2C1CCC3C4CCC(C)(C)CC4CCC23",
        "Витамин A (ретинол)": "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C",
        "Витамин B1 (тиамин)": "CC1=C(SC=N1)CN(C)C2=CN=C(N)N=C2N",
        "Витамин B2 (рибофлавин)": "CC1=CC2=C(C=C1N3C4=C(N(C=N3)C)C(=O)NC(=O)C4=C)C(C(C(C(CO)O)O)O)O",
        "Витамин B3 (ниацин)": "OC(=O)C1=CC=CN=C1",
        "Витамин B5 (пантотеновая кислота)": "CC(C)(CO)C(O)C(=O)NCCC(O)=O",
        "Витамин B6 (пиридоксин)": "CC1=C(CO)C(CO)=C(O)N=C1",
        "Витамин B7 (биотин)": "CN1C(=O)NC2C1SCC2CCC(=O)O",
        "Витамин B9 (фолиевая кислота)": "OC(=O)C(N)CCC(=O)NC1=CC=C(C2=CN=C3C(N)=NC=NC3=C2)C=C1",
        "Витамин B12 (кобаламин)": "CC1=CC2=C(C=C1N3C4=C(N(C=N3)C)C(=O)NC(=O)C4=C)C(C(C(C(CO)O)O)O)O",
        "Витамин C (аскорбиновая кислота)": "C(C(C1C(=O)C(O)=C(O)O1)O)O",
        "Витамин D2 (эргокальциферол)": "CC(C)C=CC(C)C1CCC2C1=CC=C3C2(C)CCC4C3(C)CCCC4(C)O",
        "Витамин D3 (холекальциферол)": "CC(C)CCCC(C)C1CCC2C1=CC=C3C2(C)CCC4C3(C)CCCC4(C)O",
        "Витамин E (альфа-токоферол)": "CC(C)CCCC(C)CCCC(C)CCCC1(C)CCC2CC(C)(C)CCC2(C)O1",
        "Витамин K1 (филлохинон)": "CC1=C(C)C(=O)c2c(O)c3C[C@@](C)(CCC=C(C)CCC=C(C)CCC=C(C)C)CCc3cc2C1=O",
        "Витамин K2 (менахинон)": "CC1=C(C)C(=O)c2c(O)c3C[C@@](C)(CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)CCc3cc2C1=O",
        "Бензол": "c1ccccc1",
        "Толуол": "Cc1ccccc1",
        "Ксилол": "Cc1cccc(C)c1",
        "Нафталин": "c1ccc2ccccc2c1",
        "Антрацен": "c1ccc2cc3ccccc3cc2c1",
        "Фенол": "Oc1ccccc1",
        "Анилин": "Nc1ccccc1",
        "Нитробензол": "[O-][N+](=O)c1ccccc1",
        "Хлорбензол": "Clc1ccccc1",
        "Бромбензол": "Brc1ccccc1",
        "Стирол": "C=Cc1ccccc1",
        "Бензальдегид": "O=Cc1ccccc1",
        "Бензойная кислота": "OC(=O)c1ccccc1",
        "Ацетон": "CC(C)=O",
        "Этанол": "CCO",
        "Метанол": "CO",
        "Изопропанол": "CC(C)O",
        "Бутанол": "CCCCO",
        "Диэтиловый эфир": "CCOCC",
        "Тетрагидрофуран": "C1CCOC1",
        "Диметилформамид": "CN(C)C=O",
        "Диметилсульфоксид": "CS(C)=O",
        "Ацетонитрил": "CC#N",
        "Уксусная кислота": "CC(=O)O",
        "Муравьиная кислота": "C(=O)O",
        "Пропионовая кислота": "CCC(=O)O",
        "Масляная кислота": "CCCC(=O)O",
        "Валериановая кислота": "CCCCC(=O)O",
        "Капроновая кислота": "CCCCCC(=O)O",
        "Этиленгликоль": "OCCO",
        "Глицерин": "OCC(O)CO",
        "Формальдегид": "C=O",
        "Ацетальдегид": "CC=O",
        "Пропионовый альдегид": "CCC=O",
        "Бутиловый альдегид": "CCCC=O",
        "Ацетамид": "CC(N)=O",
        "Мочевина": "NC(=O)N",
        "Гуанидин": "NC(=N)N",
        "Имидазол": "C1=CNC=N1",
        "Пиридин": "C1=CC=NC=C1",
        "Пиримидин": "C1=CC=NC=N1",
        "Пури": "C1=NC2=C(N1)N=CN2",
        "Пиррол": "C1=CC=CN1",
        "Фуран": "C1=CC=CO1",
        "Тиофен": "C1=CC=CS1",
        "Аденин": "NC1=NC=NC2=C1N=CN2",
        "Гуанин": "NC1=NC2=C(N1O)N=CN2",
        "Цитозин": "NC1=CC(=O)NC(=O)N1",
        "Тимин": "CC1=CC(=O)NC(=O)N1",
        "Урацил": "OC1=CC(=O)NC(=O)N1",
        "Аденозин": "NC1=NC=NC2=C1N=CN2[C@H]3[C@@H]([C@@H](O3)CO)O",
        "Гуанозин": "NC1=NC2=C(N1O)N=CN2[C@H]3[C@@H]([C@@H](O3)CO)O",
        "Цитидин": "NC1=CC(=O)NC(=O)N1[C@H]2[C@@H]([C@@H](O2)CO)O",
        "Уридин": "OC1=CC(=O)NC(=O)N1[C@H]2[C@@H]([C@@H](O2)CO)O",
        "Тимидин": "CC1=CC(=O)NC(=O)N1[C@H]2[C@@H]([C@@H](O2)CO)O",
        "АТФ": "NC1=NC=NC2=C1N=CN2[C@H]3[C@@H]([C@@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O",
        "АДФ": "NC1=NC=NC2=C1N=CN2[C@H]3[C@@H]([C@@H](O3)COP(=O)(O)OP(=O)(O)O)O",
        "АМФ": "NC1=NC=NC2=C1N=CN2[C@H]3[C@@H]([C@@H](O3)COP(=O)(O)O)O",
        "цАМФ": "NC1=NC=NC2=C1N=CN2[C@H]3[C@@H]([C@@H](O3)COP(=O)(O)O)O",
        "НАД": "NC(=O)c1ccc[n+](c1)[C@H]2[C@@H]([C@@H](O2)COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@H]([C@@H](O3)N4C=NC5=C4N=CN5)CO)O",
        "НАДФ": "NC(=O)c1ccc[n+](c1)[C@H]2[C@@H]([C@@H](O2)COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@H]([C@@H](O3)N4C=NC5=C4N=CN5)COP(=O)(O)O)O",
        "ФАД": "CC1=C(C)C(C)=C(C=C1N2C3=C(N(C=N3)C)C(=O)NC(=O)C2=C)C(C(C(C(CO)O)O)O)O",
        "Ресвератрол": "OC1=CC=C(C=C1)/C=C/C2=CC(O)=C(O)C=C2",
        "Куркумин": "COc1ccc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)cc1O",
        "Капсаицин": "COc1ccc(C=CC(=O)NCC(C)C)cc1O",
        "Пиперин": "O=C(NC1CCCC1)C=CC=CC2=CC=C(OC)C=C2OC",
        "Гингерол": "CCCCC(=O)C(CO)c1ccc(O)c(OC)c1",
        "Эвгенол": "COc1ccc(CC=C)cc1O",
        "Анетол": "COc1ccc(CC=C)cc1",
        "Циннамальдегид": "O=CC=Cc1ccccc1",
        "Ванилин": "COc1ccc(C=O)cc1O",
        "Салициловая кислота": "OC(=O)c1ccccc1O",
        "Бензойная кислота": "OC(=O)c1ccccc1",
        "Сорбиновая кислота": "CC=CC=CC(=O)O",
        "Бензилбензоат": "O=C(Oc1ccccc1)c2ccccc2",
        "Метилсалицилат": "COC(=O)c1ccccc1O",
        "Этилсалицилат": "CCOC(=O)c1ccccc1O",
        "Этилен": "C=C",
        "Пропилен": "CC=C",
        "Стирол": "C=Cc1ccccc1",
        "Винилхлорид": "C=CCl",
        "Акрилонитрил": "C=CC#N",
        "Метилметакрилат": "COC(=O)C(=C)C",
        "Акриловая кислота": "C=CC(=O)O",
        "Метакриловая кислота": "CC(=C)C(=O)O",
        "Винилацетат": "CC(=O)OC=C",
        "Тетрафторэтилен": "FC(F)=C(F)F",
        "Капролактам": "O=C1CCCCCN1",
        "Этиленоксид": "C1CO1",
        "Пропиленоксид": "CC1CO1",
        "Бисфенол A": "CC(C)(c1ccc(O)cc1)c2ccc(O)cc2",
        "Фталевый ангидрид": "O=C1OC(=O)c2ccccc12",
        "Терефталевая кислота": "OC(=O)c1ccc(cc1)C(=O)O",
        "Изофталевая кислота": "OC(=O)c1cccc(c1)C(=O)O",
        "Адипиновая кислота": "OC(=O)CCCCC(=O)O",
        "Себациновая кислота": "OC(=O)CCCCCCCCC(=O)O",
        "Гексаметилендиамин": "NCCCCCCN",
        "Этилендиамин": "NCCN",
        "Диэтилентриамин": "NCCNCCN",
        "Лавандулол": "CC(C)=CCC(C)(O)C=C",
        "Цитронеллол": "CC(C)CCCC(C)CCO",
        "Гераниол": "CC(C)=CCCC(C)=CCO",
        "Нерол": "CC(C)=CCCC(C)=CCO",
        "Линалоол": "CC(C)=CCCC(C)(O)C=C",
        "Цитраль": "CC(C)=CCCC(C)=CC=O",
        "Ментон": "CC(C)C1CCC(C)C(=O)C1",
        "Карвон": "CC1CCC(C)=CC1=O",
        "Пулегон": "CC1CCC(C)=CC1(=O)C",
        "Туйон": "CC1CC2CC1C(=O)C2(C)C",
        "Камфора": "CC1(C)C2CCC1(C)C(=O)C2",
        "Борнеол": "CC1(C)C2CCC1(C)C(O)C2",
        "Изокамфора": "CC1(C)C2CCC1(C)C(=O)C2",
        "Сандалол": "CC1CC2CC(C1)C2(C)CO",
        "Ветиверол": "CCC1CC2C(C1)C3=C(C2)C(=O)C=C3O",
        "Пачулиол": "CC1CC2CC(C1)C2(C)C3=CC=C(O)C=C3",
        "Жасмон": "CC1CCC(C)=CC1=O",
        "Ионон": "CC1CCC(C)=CC1=O",
        "Дамаскон": "CC1CCC(C)=CC1=O",
        "Амброксид": "CC1CC2CC(C1)C2(C)C3=CC=C(O)C=C3",
        "Стрихнин": "CN1CCC2(CC1)c3c[nH]c4cccc3C2N5C(=O)C6=CC=CC=C6C5",
        "Бруцин": "COC1=CC2=C(C=C1)C3CC4(C2)C(C3)N5CCC(C5)C=C6C7=CN=C6C=C7",
        "Кониин": "CCCC1CCCN1",
        "Никотин": "CN1CCC[C@H]1c2cccnc2",
        "Атропин": "CN1C2CCC1CC(C2)OC(=O)C(CO)C3=CC=CC=C3",
        "Скополамин": "CN1C2CCC1CC(C2)OC(=O)C(CO)C3=CC=CC=C3",
        "Кураре": "CN1CCC2(CC1)c3c[nH]c4cccc3C2N5C(=O)C6=CC=CC=C6C5",
        "Тетродотоксин": "OC1C(O)C2(O)C3C4C(O)C5(O)C6(O)C(O)C7(O)C(O)C(O)C(O)C(O)C12C34C56C7",
        "Рицин": "NC(=O)C1C(O)C(O)C(O)C(O)C1O",
        "Афлатоксин": "COC1=C2C(=O)C3=C(C=C(C)C3)OC2=CC4=C1C5=C(C=C4)OC(=O)C5",
        "Ботулотоксин": "CC(C)C(N)C(=O)N",
        "Диоксин": "Clc1c(Cl)cc2c(c1Cl)oc3c(Cl)c(Cl)cc(Cl)c23",
        "Зарин": "CC(C)OP(=O)(F)OC",
        "Зоман": "CC(C)C(C)OP(=O)(F)OC",
        "VX": "CCOP(=O)(OCC)SCCN(C)C",
        "Фосген": "ClC(=O)Cl",
        "Иприт": "ClCCSCCCl",
        "Хлорпикрин": "ClC(Cl)(Cl)[N+](=O)[O-]",
        "Люизит": "ClC=C(Cl)AsCl2",
        "ДДТ": "Clc1ccc(C(c2ccc(Cl)cc2)(c2ccc(Cl)cc2)Cl)cc1",
        "Диэльдрин": "ClC1C(Cl)C2C3C4C1C5C4C(Cl)C(Cl)C5C3(Cl)O2",
        "Эндрин": "ClC1C(Cl)C2C3C4C1C5C4C(Cl)C(Cl)C5C3(Cl)O2",
        "Хлордан": "ClC1C(Cl)C2C3C4C1C5C4C(Cl)C(Cl)C5C3(Cl)O2",
        "Токсафен": "ClC1C(Cl)C2C3C4C1C5C4C(Cl)C(Cl)C5C3(Cl)O2",
        "ПХБ": "Clc1c(Cl)cc2c(c1Cl)cc3c(Cl)c(Cl)cc(Cl)c23",
        "Бензпирен": "c1ccc2cc3cc4ccc5cc6ccc7ccccc7c6c5c4cc3cc2c1",
        "Акриламид": "C=CC(N)=O",
        "Бисфенол A": "CC(C)(c1ccc(O)cc1)c2ccc(O)cc2",
        "Фталаты": "CCOC(=O)c1ccccc1C(=O)OCC",
        "Перхлорэтилен": "ClC(Cl)=C(Cl)Cl",
        "Трихлорэтилен": "ClC(Cl)=CCl",
        "Винилхлорид": "C=CCl",
        "Бензол": "c1ccccc1",
        "Толуол": "Cc1ccccc1",
        "Ксилол": "Cc1cccc(C)c1",
        "Нафталин": "c1ccc2ccccc2c1",
        "Антрацен": "c1ccc2cc3ccccc3cc2c1",
        "Фенантрен": "c1ccc2cc3ccccc3cc2c1",
        "Трифенилфосфин": "P(c1ccccc1)(c2ccccc2)c3ccccc3",
        "Триэтиламин": "CCN(CC)CC",
        "Пиридин": "C1=CC=NC=C1",
        "Имидазол": "C1=CNC=N1",
        "Триэтилсилилхлорид": "CC[Si](Cl)(CC)CC",
        "Трет-бутилдиметилсилилхлорид": "CC[Si](Cl)(C)C(C)(C)C",
        "Диизопропилазодикарбоксилат": "CC(C)OC(=O)N=N C(=O)OC(C)C",
        "Трифенилметилхлорид": "ClC(c1ccccc1)(c2ccccc2)c3ccccc3",
        "Бензилбромид": "Brc1ccccc1",
        "Аллилбромид": "BrCC=C",
        "Пропаргилбромид": "BrCC#C",
        "Метилйодид": "CI",
        "Этилиодид": "CCI",
        "Бензилйодид": "Ic1ccccc1",
        "Триметилсилилйодид": "C[Si](I)(C)C",
        "Тетрабутиламмонийбромид": "CCCC[N+](CCCC)(CCCC)CCCC.[Br-]",
        "Тетрабутиламмониййодид": "CCCC[N+](CCCC)(CCCC)CCCC.[I-]",
        "Тетраметиламмонийхлорид": "C[N+](C)(C)C.[Cl-]",
        "Тетраэтиламмонийхлорид": "CC[N+](CC)(CC)CC.[Cl-]",
        "Палладий на угле": "[Pd]",
        "Платина на угле": "[Pt]",
        "Никель Ренея": "[Ni]",
        "Тетракис(трифенилфосфин)палладий": "P(c1ccccc1)(c2ccccc2)(c3ccccc3)[Pd]P(c4ccccc4)(c5ccccc5)c6ccccc6",
        "Дихлорбис(трифенилфосфин)палладий": "Cl[Pd](Cl)P(c1ccccc1)(c2ccccc2)c3ccccc3",
        "Ацетат палладия": "CC(=O)O[Pd]OC(=O)C",
        "Хлорид палладия": "Cl[Pd]Cl",
        "Бромид палладия": "Br[Pd]Br",
        "Триэтиламин": "CCN(CC)CC",
        "Диизопропилэтиламин": "CCN(CC(C)C)C(C)C",
        "1,8-Диазабицикло[5.4.0]ундец-7-ен": "C1CN2CCC1CC2",
        "4-Диметиламинопиридин": "CN(C)c1ccnc(c1)N",
        "Имидазол": "C1=CNC=N1",
        "1-Метилимидазол": "Cn1ccnc1",
        "Триэтилборан": "CC[B](CC)CC",
        "Триэтилалюминий": "CC[Al](CC)CC",
        "Диэтилцинк": "CC[Zn]CC",
        "Триметилалюминий": "C[Al](C)C",
        "Триметилборан": "C[B](C)C",
        "Вода": "O",
        "Метанол": "CO",
        "Этанол": "CCO",
        "Изопропанол": "CC(C)O",
        "н-Бутанол": "CCCCO",
        "трет-Бутанол": "CC(C)(C)O",
        "Этиленгликоль": "OCCO",
        "Глицерин": "OCC(O)CO",
        "Диметилсульфоксид": "CS(C)=O",
        "Диметилформамид": "CN(C)C=O",
        "Ацетонитрил": "CC#N",
        "Ацетон": "CC(C)=O",
        "Метилэтилкетон": "CCC(C)=O",
        "Циклогексанон": "O=C1CCCCC1",
        "Тетрагидрофуран": "C1CCOC1",
        "Диоксан": "C1COCCO1",
        "Диэтиловый эфир": "CCOCC",
        "Трет-бутилметиловый эфир": "CC(C)(C)OCC",
        "Дихлорметан": "ClCCl",
        "Хлороформ": "ClC(Cl)Cl",
        "Четыреххлористый углерод": "ClC(Cl)(Cl)Cl",
        "1,2-Дихлорэтан": "ClCCl",
        "Трихлорэтилен": "ClC(Cl)=CCl",
        "Перхлорэтилен": "ClC(Cl)=C(Cl)Cl",
        "Бензол": "c1ccccc1",
        "Толуол": "Cc1ccccc1",
        "Ксилол": "Cc1cccc(C)c1",
        "Гексан": "CCCCCC",
        "Гептан": "CCCCCCC",
        "Октан": "CCCCCCCC",
        "Циклогексан": "C1CCCCC1",
        "Циклопентан": "C1CCCC1",
        "Пиридин": "C1=CC=NC=C1",
        "Уксусная кислота": "CC(=O)O",
        "Муравьиная кислота": "C(=O)O",
        "Триэтиламин": "CCN(CC)CC",
        "Диэтиламин": "CCNCC",
        "Трифторметансульфокислота": "OS(=O)(=O)C(F)(F)F",
        "Трифторметансульфонилтрифлат": "FC(F)(F)S(=O)(=O)OS(=O)(=O)C(F)(F)F",
    }

  
    
    selected_name = st.selectbox("Select Molecule:", list(examples.keys()))
    smiles_input = st.text_input("SMILES Code:", value=examples[selected_name])
    analyze_clicked = st.button("❤️ Analyze Structure")
    
    st.markdown("---")
    st.caption("Powered by RDKit & Streamlit")

if analyze_clicked:
    with st.spinner('Processing data...'):
        data = calculate_properties(smiles_input)
    
    if "error" in data:
        st.error(f"🚫 Invalid Structure: {data['error']}")
        st.stop()
    
    mol = data.pop("mol")
    violations = data.pop("Violations")
    inchi_key = data.pop("InChIKey")
    
    c1, c2, c3, c4 = st.columns(4)

# Красный цвет для метрик (HTML + CSS)
metric_style = """
<style>
.red-metric {
    background: white;
    border-radius: 12px;
    padding: 16px;
    text-align: center;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    border: 1px solid #efefef;
}
.red-metric-label {
    color: #e63946 !important;
    font-size: 0.9rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 8px;
}
.red-metric-value {
    color: #e63946 !important;
    font-size: 2rem;
    font-weight: 700;
}
</style>
"""

st.markdown(metric_style, unsafe_allow_html=True)

with c1:
    st.markdown(
        f'<div class="red-metric"><div class="red-metric-label">MW</div><div class="red-metric-value">{data["MW"]}</div></div>',
        unsafe_allow_html=True
    )
with c2:
    st.markdown(
        f'<div class="red-metric"><div class="red-metric-label">LogP</div><div class="red-metric-value">{data["LogP"]}</div></div>',
        unsafe_allow_html=True
    )
with c3:
    st.markdown(
        f'<div class="red-metric"><div class="red-metric-label">TPSA</div><div class="red-metric-value">{data["TPSA"]}</div></div>',
        unsafe_allow_html=True
    )
with c4:
    lipinski_value = "✅ PASS" if violations == 0 else f"⚠️ {violations}"
    st.markdown(
        f'<div class="red-metric"><div class="red-metric-label">Lipinski</div><div class="red-metric-value">{lipinski_value}</div></div>',
        unsafe_allow_html=True
    )

    st.markdown("---")
    
    tab_struct, tab_3d, tab_data = st.tabs(["📸 Structure", "🧊 3D View", "📝 Details"])
    
    with tab_struct:
        col_img, col_info = st.columns([1, 1])
        with col_img:
            img = MolToImage(mol, size=(400, 400))
            st.image(img, width=400)
            
        with col_info:
            st.subheader("Chemical Identity")
            st.code(f"SMILES: {smiles_input}", language="bash")
            st.code(f"InChIKey: {inchi_key}", language="bash")
            
            st.markdown("### 💡 Insight")
            if data['LogP'] > 3:
                st.info("🌊 Hydrophobic molecule. Good for membrane permeability.")
            else:
                st.info("💧 Hydrophilic tendencies. Good solubility expected.")
                
            if data['MW'] < 300:
                st.success("🪶 Lightweight fragment. Good for drug optimization.")

    with tab_3d:
        st.markdown("#### Interactive Model")
        html_3d = generate_3d_html(mol, width=700, height=450)
        components.html(html_3d, height=470, scrolling=False)
        
    with tab_data:
        df = pd.DataFrame({
            "Property": ["H-Bond Donors", "H-Bond Acceptors", "Rotatable Bonds", "Rule of 5 Violations"],
            "Value": [data['HBD'], data['HBA'], data['RotBonds'], violations]
        })
        st.table(df)
        
        st.markdown("### 📥 Download Assets")
        col_d1, col_d2 = st.columns(2)
        
        buf = BytesIO()
        img.save(buf, format='PNG')
        with col_d1:
            st.download_button(
                label="Download Image (PNG)",
                data=buf.getvalue(),
                file_name=f"{inchi_key[:8]}.png",
                mime="image/png"
            )
        
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        sdf = Chem.MolToMolBlock(mol_3d)
        with col_d2:
            st.download_button(
                label="Download Structure (SDF)",
                data=sdf,
                file_name=f"{inchi_key[:8]}.sdf",
                mime="chemical/x-mdl-sdfile"
            )

else:
    st.markdown("""
    <div style='display:flex; justify-content:center; align-items:center; height:400px; flex-direction:column; color:#dbdbdb;'>
        <h1 style='font-size:4rem; margin:0;'>🧪</h1>
        <h3 style='color:#8e8e8e;'>Ready to discover?</h3>
        <p>Select a molecule in the sidebar to start.</p>
    </div>
    """, unsafe_allow_html=True)
