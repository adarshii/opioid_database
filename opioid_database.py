import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import pandas as pd
import requests
import pubchempy as pcp
import base64
from io import BytesIO

# Sample opioid drug database
OPIOID_DATA = [
    {
        "name": "Morphine",
        "iupac": "(5α,6α)-7,8-Didehydro-4,5-epoxy-17-methylmorphinan-3,6-diol",
        "smiles": "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        "formula": "C17H19NO3",
        "weight": 285.34,
        "category": "Natural Opioid",
        "uses": ["Severe pain management", "Palliative care"],
        "toxicity": "High - Respiratory depression, addiction potential",
        "symptoms": ["Drowsiness", "Constipation", "Nausea", "Respiratory depression", "Coma"],
        "side_effects": ["Dependence", "Tolerance", "Withdrawal symptoms"],
        "precautions": ["Avoid alcohol", "Monitor respiratory function", "Risk of addiction"],
        "pubchem_cid": "5288826"
    },
    {
        "name": "Fentanyl",
        "iupac": "N-(1-(2-Phenylethyl)-4-piperidinyl)-N-phenylpropanamide",
        "smiles": "CCN(CC1=CC=CC=C1)C(=O)C2CCN(CC2)CCC3=CC=CC=C3",
        "formula": "C22H28N2O",
        "weight": 336.47,
        "category": "Synthetic Opioid",
        "uses": ["Anesthesia", "Chronic pain management"],
        "toxicity": "Extremely High - 50-100x more potent than morphine",
        "symptoms": ["Respiratory arrest", "Pinpoint pupils", "Unconsciousness"],
        "side_effects": ["Muscle rigidity", "Hypotension", "Bradycardia"],
        "precautions": ["For medical professionals only", "Risk of fatal overdose"],
        "pubchem_cid": "3345"
    },
    {
        "name": "Oxycodone",
        "iupac": "(5R,9R,13S,14S)-4,5α-Epoxy-14-hydroxy-3-methoxy-17-methylmorphinan-6-one",
        "smiles": "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(CC4=O)O",
        "formula": "C18H21NO4",
        "weight": 315.36,
        "category": "Semi-Synthetic Opioid",
        "uses": ["Moderate to severe pain", "Post-surgical pain"],
        "toxicity": "High - Abuse potential, respiratory depression",
        "symptoms": ["Dizziness", "Confusion", "Hypotension", "Respiratory depression"],
        "side_effects": ["Constipation", "Nausea", "Sedation", "Addiction"],
        "precautions": ["Avoid with alcohol", "Monitor for dependence", "Risk of misuse"],
        "pubchem_cid": "5284603"
    },
    {
        "name": "Hydrocodone",
        "iupac": "(4R,4aR,7S,7aR,12bS)-3-Methoxy-9-hydroxy-2,4,4a,7,7a,13-hexahydro-1H-4,12-methanobenzofuro[3,2-e]isoquinoline-7-one",
        "smiles": "CN1CCC23C4C1C5=C2C(=C(C=C5)OC)OC3C(CC4=O)O",
        "formula": "C18H21NO3",
        "weight": 299.37,
        "category": "Semi-Synthetic Opioid",
        "uses": ["Pain relief", "Cough suppression"],
        "toxicity": "High - Similar to oxycodone",
        "symptoms": ["Sedation", "Respiratory depression", "Dizziness"],
        "side_effects": ["Constipation", "Nausea", "Headache", "Dependence"],
        "precautions": ["Avoid driving", "Risk of addiction", "Monitor for misuse"],
        "pubchem_cid": "5284569"
    },
    {
        "name": "Codeine",
        "iupac": "(5α,6α)-7,8-Didehydro-4,5-epoxy-3-methoxy-17-methylmorphinan-6-ol",
        "smiles": "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O",
        "formula": "C18H21NO3",
        "weight": 299.36,
        "category": "Natural Opioid",
        "uses": ["Mild to moderate pain", "Cough suppression"],
        "toxicity": "Moderate - Less potent than morphine",
        "symptoms": ["Drowsiness", "Constipation", "Nausea"],
        "side_effects": ["Dependence with long-term use", "Allergic reactions"],
        "precautions": ["Avoid with CNS depressants", "Risk of respiratory depression"],
        "pubchem_cid": "5284371"
    },
    {
        "name": "Methadone",
        "iupac": "(6-Dimethylamino-4,4-diphenylheptan-3-one)",
        "smiles": "CCC(=O)C(CC(C)N(C)C)(C1=CC=CC=C1)C2=CC=CC=C2",
        "formula": "C21H27NO",
        "weight": 309.45,
        "category": "Synthetic Opioid",
        "uses": ["Opioid addiction treatment", "Chronic pain management"],
        "toxicity": "High - Long half-life, risk of accumulation",
        "symptoms": ["QT prolongation", "Respiratory depression", "Cardiac arrhythmias"],
        "side_effects": ["Constipation", "Sweating", "Sexual dysfunction"],
        "precautions": ["Requires medical supervision", "Risk of overdose", "Monitor ECG"],
        "pubchem_cid": "4095"
    }
]

# Utility functions
def get_pubchem_link(cid):
    return f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

def get_compound_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(400, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode()
    return None

def get_3d_structure(cid):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d&response_type=save"
        response = requests.get(url)
        return response.content
    except:
        return None

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {}
    
    return {
        "Molecular Weight": f"{Descriptors.MolWt(mol):.2f} g/mol",
        "LogP": f"{Descriptors.MolLogP(mol):.2f}",
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Polar Surface Area": f"{Descriptors.TPSA(mol):.2f} Å²"
    }

# Home page
def home_page():
    st.title("Opioid Drugs Database")
    st.image("https://images.unsplash.com/photo-1614064641938-3bbee52942c7?auto=format&fit=crop&w=1200", 
             use_column_width=True, caption="Opioid Molecular Structure")
    
    st.markdown("""
    ## Understanding Opioid Drugs
    
    Opioids are a class of drugs that include:
    - Natural opioid analgesics (e.g., morphine, codeine)
    - Semi-synthetic opioids (e.g., oxycodone, hydrocodone)
    - Fully synthetic opioids (e.g., fentanyl, methadone)
    
    ### Medical Uses:
    - Management of acute and chronic pain
    - Anesthesia adjuncts
    - Cough suppression
    - Treatment of opioid use disorder
    
    ### Risks and Concerns:
    - High potential for addiction and dependence
    - Respiratory depression (can be fatal)
    - Increasing rates of overdose deaths
    - Diversion and misuse
    
    **Use this database to access scientific information about opioid medications.**
    """)
    
    st.warning("**Important:** This information is for educational purposes only. Never use prescription medications without medical supervision.")
    
    # Search functionality
    st.subheader("Search Opioid Database")
    drug_names = [drug["name"] for drug in OPIOID_DATA]
    drug_names.sort()
    selected_drug = st.selectbox("Select or type an opioid drug:", drug_names)
    
    if st.button("View Drug Details"):
        st.session_state.current_drug = selected_drug
        st.experimental_rerun()

# Drug detail page
def drug_detail_page():
    drug_name = st.session_state.get("current_drug", "")
    
    if not drug_name:
        st.error("No drug selected. Please go back to home page.")
        return
    
    # Find drug in database
    drug = next((item for item in OPIOID_DATA if item["name"] == drug_name), None)
    
    if not drug:
        st.error("Drug information not available")
        return
    
    st.title(f"{drug_name} - Opioid Drug Profile")
    st.markdown(f"**Category:** {drug['category']}")
    
    # Create tabs for different sections
    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["Overview", "Chemical Structure", "Pharmacology", "Safety", "References"]
    )
    
    with tab1:  # Overview tab
        col1, col2 = st.columns([1, 2])
        
        with col1:
            # Display 2D structure
            img_data = get_compound_image(drug["smiles"])
            if img_data:
                st.image(f"data:image/png;base64,{img_data}", 
                         caption=f"2D Structure of {drug_name}")
            else:
                st.warning("Structure image not available")
            
            # Basic properties
            st.subheader("Molecular Properties")
            st.markdown(f"**IUPAC Name:** {drug['iupac']}")
            st.markdown(f"**Molecular Formula:** {drug['formula']}")
            st.markdown(f"**Molecular Weight:** {drug['weight']} g/mol")
            
            # PubChem link
            st.subheader("External References")
            st.markdown(f"[PubChem Entry]({get_pubchem_link(drug['pubchem_cid'])})")
        
        with col2:
            # Medical uses
            st.subheader("Medical Uses")
            for use in drug["uses"]:
                st.markdown(f"- {use}")
            
            # Precautions
            st.subheader("Precautions")
            for precaution in drug["precautions"]:
                st.markdown(f"⚠️ {precaution}")
    
    with tab2:  # Chemical Structure tab
        st.header("Chemical Structure Information")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Display 2D structure
            st.subheader("2D Molecular Structure")
            img_data = get_compound_image(drug["smiles"])
            if img_data:
                st.image(f"data:image/png;base64,{img_data}", 
                         use_column_width=True,
                         caption=f"2D Structure of {drug_name}")
            else:
                st.warning("2D structure image not available")
            
            # SMILES notation
            st.subheader("SMILES Notation")
            st.code(drug["smiles"], language="text")
            
            # Calculated properties
            st.subheader("Physicochemical Properties")
            properties = calculate_properties(drug["smiles"])
            for prop, value in properties.items():
                st.markdown(f"**{prop}:** {value}")
        
        with col2:
            # 3D structure placeholder
            st.subheader("3D Molecular Structure")
            st.markdown(f"""
            <div style="border:1px solid #ccc; padding:20px; border-radius:10px; background-color:#f9f9f9; height:400px">
                <h4 style="color:#555">Interactive 3D Structure</h4>
                <p>For full 3D visualization, visit the PubChem entry:</p>
                <a href="{get_pubchem_link(drug['pubchem_cid'])}" target="_blank">
                    <button style="background-color:#1e3a8a; color:white; border:none; padding:10px 20px; border-radius:5px; cursor:pointer">
                        View on PubChem
                    </button>
                </a>
                <p style="margin-top:20px">PubChem CID: {drug['pubchem_cid']}</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Download SDF file
            sdf_data = get_3d_structure(drug["pubchem_cid"])
            if sdf_data:
                st.download_button(
                    label="Download 3D Structure (SDF)",
                    data=sdf_data,
                    file_name=f"{drug_name}_3d_structure.sdf",
                    mime="chemical/x-mdl-sdfile"
                )
            else:
                st.warning("3D structure data not available for download")
    
    with tab3:  # Pharmacology tab
        st.header("Pharmacology and Clinical Use")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Therapeutic Uses")
            for use in drug["uses"]:
                st.markdown(f"- {use}")
            
            st.subheader("Mechanism of Action")
            st.markdown("""
            Opioids work by binding to specific receptors in the brain and spinal cord:
            - μ-opioid receptors (MOR) - primary analgesic effect
            - δ-opioid receptors (DOR)
            - κ-opioid receptors (KOR)
            
            This binding:
            - Inhibits neurotransmitter release
            - Hyperpolarizes neurons
            - Modulates pain transmission
            """)
            
        with col2:
            st.subheader("Pharmacokinetics")
            st.markdown("""
            **Typical Pharmacokinetic Properties:**
            - Onset of action: 5-30 minutes (IV) to 30-90 minutes (oral)
            - Duration of effect: 3-6 hours for most opioids
            - Metabolism: Primarily hepatic (CYP450 system)
            - Excretion: Renal
            
            **Note:** Actual parameters vary significantly between specific opioids
            """)
            
            st.subheader("Dosage Information")
            st.warning("Dosage must be determined by a qualified medical professional based on:")
            st.markdown("""
            - Patient's pain level
            - Previous opioid exposure
            - Medical condition
            - Other medications
            - Risk factors for respiratory depression
            """)
    
    with tab4:  # Safety tab
        st.header("Safety and Toxicity Profile")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Toxicity Information")
            st.error(f"**Toxicity Level:** {drug['toxicity']}")
            
            st.subheader("Overdose Symptoms")
            for symptom in drug["symptoms"]:
                st.markdown(f"- {symptom}")
            
            st.subheader("Management of Overdose")
            st.markdown("""
            1. Call emergency services immediately
            2. Administer naloxone if available
            3. Provide respiratory support
            4. Monitor vital signs continuously
            """)
            
        with col2:
            st.subheader("Side Effects")
            st.warning("Common side effects include:")
            for effect in drug["side_effects"]:
                st.markdown(f"- {effect}")
            
            st.subheader("Contraindications")
            st.markdown("""
            - Significant respiratory depression
            - Acute or severe bronchial asthma
            - Known or suspected gastrointestinal obstruction
            - Hypersensitivity to specific opioid
            """)
            
            st.subheader("Drug Interactions")
            st.markdown("""
            - **CNS depressants:** Alcohol, benzodiazepines (increased sedation)
            - **MAO inhibitors:** Risk of serotonin syndrome
            - **Anticholinergics:** Increased risk of constipation
            """)
    
    with tab5:  # References tab
        st.header("References and Further Reading")
        
        st.subheader("Primary Reference")
        st.markdown(f"""
        - **PubChem Entry:** [{drug_name} (CID: {drug['pubchem_cid']})]({get_pubchem_link(drug['pubchem_cid'])})
        """)
        
        st.subheader("Clinical Guidelines")
        st.markdown("""
        - [CDC Guideline for Prescribing Opioids for Chronic Pain](https://www.cdc.gov/opioids/providers/prescribing/guideline.html)
        - [WHO Guidelines for the Pharmacological Treatment of Persisting Pain in Children](https://www.who.int/publications/i/item/9789241548120)
        - [ASAM National Practice Guideline for the Treatment of Opioid Use Disorder](https://www.asam.org/quality-care/clinical-guidelines/national-practice-guideline)
        """)
        
        st.subheader("Scientific Literature")
        st.markdown("""
        - [Opioid Pharmacology Review - Nature Reviews Drug Discovery](https://www.nature.com/articles/nrd3438)
        - [Opioid Overdose Crisis - New England Journal of Medicine](https://www.nejm.org/doi/full/10.1056/nejmra1508491)
        - [Molecular Mechanisms of Opioid Action - Pharmacological Reviews](https://pharmrev.aspetjournals.org/content/73/4/1638)
        """)
        
        st.subheader("Additional Resources")
        st.markdown("""
        - [National Institute on Drug Abuse (NIDA) - Opioids](https://www.drugabuse.gov/drug-topics/opioids)
        - [SAMHSA Opioid Overdose Prevention Toolkit](https://store.samhsa.gov/product/Opioid-Overdose-Prevention-Toolkit/SMA18-4742)
        - [FDA Opioid Medications Information](https://www.fda.gov/drugs/information-drug-class/opioid-medications)
        """)
    
    # Back button
    st.divider()
    if st.button("Back to Home Page"):
        del st.session_state.current_drug
        st.experimental_rerun()

# Main app function
def main():
    # Set page config
    st.set_page_config(
        page_title="Opioid Drugs Database",
        page_icon="💊",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Initialize session state
    if "current_drug" not in st.session_state:
        st.session_state.current_drug = None
    
    # Navigation
    if st.session_state.current_drug:
        drug_detail_page()
    else:
        home_page()

if __name__ == "__main__":
    main()
