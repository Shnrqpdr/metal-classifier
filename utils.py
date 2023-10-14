import ase.db # https://wiki.fysik.dtu.dk/ase/ase/db/db.html
from pymatgen.core.composition import *
import pandas as pd
import numpy as np

def create_df_stoichiometry(stoichiometry):
    target = 'ehull'
    # Initialize DataFrame columns
    columns = [
        f'Material',
        f'Space_Group',
        f'{target}',
        f'Atom1',
        f'Atom2'
    ]
    
    if stoichiometry == 'ABC':
        columns.append(f'Atom3')

    df = pd.DataFrame(columns=columns)
    
    data = ase.db.connect('c2db-2021-06-24.db')
    data_nonMag = data.select(is_magnetic=False)
    
    for row in data_nonMag:
        try:
            formula_information = Composition(row.formula).as_dict()
            list_ele = list(formula_information.items())
            
            if stoichiometry == 'ABC':
                if len(list_ele) == 3 and all(value == 1 for key, value in list_ele):
                    new_entry = [row.formula, row.spacegroup, row[target], list_ele[0][0], list_ele[1][0], list_ele[2][0]]
                    df.loc[len(df)] = new_entry
            else:
                atom_ratio = int(stoichiometry[-1])
                if len(list_ele) == 2 and list_ele[0][1] == 1 and list_ele[1][1] == atom_ratio:
                    new_entry = [row.formula, row.spacegroup, row[target], list_ele[0][0], list_ele[1][0]]
                    df.loc[len(df)] = new_entry
                    
        except Exception as e:
            print(f"Caught an exception: {e}")

    return df

def calculate_statistical_features(df_materials, df_atoms, atom_columns):
    # Initialize an empty DataFrame to hold the new features
    df_features = pd.DataFrame()
    
    # List to hold all statistical feature names
    stats_columns = ['media', 'max', 'min', 'deviation']
    
    # Loop through each atomic property in df_atoms
    for col in df_atoms.columns:
        if col == 'Element':
            continue
        # Create empty lists to hold statistical values for each material
        stats_values = {stat: [] for stat in stats_columns}
        
        # Loop through each row in df_materials
        for index, row in df_materials.iterrows():
            # Extract the atomic properties for the atoms in this material
            atom_properties = [df_atoms.loc[df_atoms['Element'] == row[atom_col], col].values[0] for atom_col in atom_columns if row[atom_col] in df_atoms['Element'].values]
            # Skip if no valid atomic properties found
            if len(atom_properties) == 0:
                continue
            
            # Calculate the statistical features for this material and property
            stats_values['media'].append(np.mean(atom_properties))
            stats_values['max'].append(np.max(atom_properties))
            stats_values['min'].append(np.min(atom_properties))
            stats_values['deviation'].append(np.std(atom_properties))
        
        # Add these statistical features to df_features
        for stat, values in stats_values.items():
            df_features[f'{stat}_{col}'] = values
            
    # Add these new features to the original df_materials DataFrame and return
    df_materials = pd.concat([df_materials.reset_index(drop=True), df_features.reset_index(drop=True)], axis=1)
    
    return df_materials