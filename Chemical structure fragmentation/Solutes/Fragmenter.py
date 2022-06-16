#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import operator
import rdkit
from fragmenter import fragmenter
from rdkit import Chem


# In[ ]:


def info_to_CSV(inchikey, SMILES, pubchem_id, fragmentation):
    
    fragmentation_array = []
    for group_number, amount in fragmentation.items():
        fragmentation_array.append(str(group_number) + ":" + str(amount))
    
    return inchikey + "," + SMILES + "," + pubchem_id + "," + "|".join(fragmentation_array)


# In[ ]:


def CSV_to_info(CSV_line, has_fragmentation = False):
    CSV_line = CSV_line.replace('\n', '')
    array = CSV_line.split(',')
    
    fragmentation = {}
    
    if has_fragmentation:
        fragmentation_array = array[3].split('|')
        for match_str in fragmentation_array:
            array2 = match_str.split(':')
            group_number = int(array2[0])
            amount = int(array2[1])
            
            fragmentation[group_number] = amount
        
    return array[0], array[1], array[2], fragmentation


# In[ ]:


def function_to_choose_fragmentation(fragmentations):
    fragmentations_descriptors = {}
    i = 0
    for fragmentation in fragmentations:
        fragmentations_descriptors[i] = [len(fragmentation)]
        i += 1
    
    sorted_fragmentations_dict = sorted(fragmentations_descriptors.items(), key=operator.itemgetter(1))

    return fragmentations[sorted_fragmentations_dict[0][0]]


# In[ ]:


def is_fragmentation_equal_to_other_fragmentation(fragmentation, other_fragmentation):
    
    for group_number, amount in fragmentation.items():
        if group_number in other_fragmentation:
            if fragmentation[group_number] != other_fragmentation[group_number]:
                return False
    return True


# In[ ]:


def log_structure_results(f, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, status = ''):
    
    f.write('https://pubchem.ncbi.nlm.nih.gov/compound/' + pubchem_id + '#section=2D-Structure\n')
    f.write(SMILES + '\n')
    f.write(inchikey + '\n')
    f.write('\n' + 'Fragmentation was successful: ' + str(success) + '\n')
    
    if status != '':
        f.write(status + '\n')
    
    if success:
        f.write('Fragmentation from the algorithm:\n')
        sorted_group_number = sorted(fragmentation.keys())
        
        for group_number in sorted_group_number:
            f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(fragmentation[group_number]).ljust(8, ' ') + '\n')
    
    f.write('\n')
    
    if len(fragmentation_reference_DB) > 0:
        f.write('Fragmentation from the reference database:\n')
        sorted_group_number = sorted(fragmentation_reference_DB.keys())
        
        for group_number in sorted_group_number:
            f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(fragmentation_reference_DB[group_number]).ljust(8, ' ') + '\n')
    
        
    f.write('\n\n')   

UNIFAC_SMARTS =  [
("CH3", "[CH3;X4]"),
("CH2", "[CH2;X4]"),
("CH", "[CH1;X4]"),
("C", "[CH0;A;X4;!R]"),
("CH2=CH", "[CH2]=[CH]"),
("CH=CH", "[CH]=[CH]"),
("CH2=C", ["[CH2]=[C]", "[CH2]=[c]"]),
("CH=C", ["[CH]=[CH0]", "[CH]=[cH0]"]),
("CH#C", "C#C"),
("PHCH", "[cH]"),
("PHC", "[cH0;!+]"),
("PHCCH3", "[c][CH3;X4]"),
("PHCCH2", "[c][CH2;X4]"),
("PHCCH", "[c][CH;X4]"),
("OH", "[OH]"),
("CH3OH", "[CH3][OH]"),
("CH3CO", "[CH3][CH0]=O"),    
("CH2CO", "[CH2][CH0]=O"),
("CHO", "[CH]=O"),
("CH2O", "[CH2]=O"),
("CH3COO", "[CH3]C(=O)[OH0;!-]"),    
("CH2COO", "[CH2]C(=O)[OH0;!-]"),
("OCH3", "[CH3][OH0]"),    
("OCH2", "[CH2][OH0]"),     
("OCH", "[CH][OH0]"),
("SI(1)3", "[C]-[Si](-[C])(-[C])-[C]"),
("OBO", "[O]-[B]-[O]"),
("SO3", "[O]=[S](-[O])=[O]"),
("CNC", "[C]-[N;X3]-[C]"),
("CN", "[C]#[N]"),
("THF", "[CH2;R][OH0]"),
("CH2N", "[CH2][N;!+]"),    
("CH3CN", "[CH3]C#N"),    
("CH2CN", "[CH2]C#N"),
("COOH", "C(=O)[OH]"),
("CH2Cl", "[CH2]Cl"),
("CH2Cl2", "[CH2](Cl)Cl"),
("CHCl3", "[CH](Cl)(Cl)Cl"),    
("CCl4", "C(Cl)(Cl)(Cl)(Cl)"),    
("CH3NO2", "[CH3][N+](=O)[O-]"),    
("CH2NO2", "[CH2][N+](=O)[O-]"),    
("thiophene", "[cH]1[cH][s;X2][cH][cH]1"),
("trifluoroethanol", "OCC(F)(F)F"),
("halothane", "FC(F)(F)C(Cl)Br"),
("C4H3S", ['[c]1[cH][s;X2][cH][cH]1', 
           '[cH]1[c][s;X2][cH][cH]1']),    
("C4H2S", ['[c]1[c][s;X2][cH][cH]1', 
           '[c]1[cH][s;X2][cH][c]1', 
           '[cH]1[c][s;X2][c][cH]1', 
           '[cH]1[c][s;X2][cH][c]1']),
("pyridine", "[n;!+]1[cH][cH][cH][cH][cH]1"),    
("C5H4N", ['[n;!+]1[c][cH][cH][cH][cH]1', 
           '[n;!+]1[cH][c][cH][cH][cH]1', 
           '[n;!+]1[cH][cH][c][cH][cH]1']),    
("H2O", "[OH2]"),
("[IM-1]", "C[NH+]1C=CN=C1"),
("[IM-1,(a)]", ['[n;X3]1([cH][n+;X3](cc1)-[CH3])', 
              '[n;X3]1([cH][n+;X3]cc1)-[CH3]']),
("[IM-(a),(b)]", "[n]1[c][n+;X3](cc1)"),
("[IM-1,(a),(b)]", ['[n;X3]1([c;X3][n+;X3](cc1)-[CH3])',
              '[n;X3]1([c;X3][n+;X3]cc1)-[CH3]',
              'n1c([n+]cc1)-[CH3]']),
("[PY-1,(a)]", ['c1ccc(c[n+;X3]1)-[CH3]', 
         'c1cccc([n+;X3]1)-[CH3]', 
         'c1cc(cc[n+;X3]1)-[CH3]']), 
("[PY]", "c1cccc[n+]1"),    
("[PYR-1,(a)]", "[N+;X4]-1(-[C]-[C]-[C]-[C]-1)(-[CH3])"),
("[PYR]", "[N+;X4]-1(-[C]-[C]-[C]-[C]-1)"),
("[PIP-1,(a)]", "[C]-1-[N+;X4](-[C]-[C]-[C]-[C]-1)(-[CH3])"),    
("[QUINUC]", "[C]-1-[N+;X4](-[C]-[C]-[C]-[C]-1)"),     
("[MO-1,(a)]", "[C]-1-[N+;X4](-[C]-[C]-[O]-[C]-1)(-[CH3])"),    
("[QUINI]", "c1ccc2c(c1)c[n+;X3](cc2)"),
("[CPP]", "[c+]1cc1"),
("[MTBDH]", "CN1CCCN2CCC[NH+]=C12"),
("[CH3N]", "[N+;!R;X4](-[CH3])"),
("[C2H5N]", "[N+;!R;X4](-[CH2]-[CH3])"),
("[C3H7N]", "[N+;!R;X4]-[CH2]-[CH2]-[CH3]"),
("[C4H9N]", "[N+;!R;X4](-[CH2]-[CH2]-[CH2]-[CH3])"),
("[C6H13N]", "[N+;!R;X4](-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3])"), 
("[C8H17N]", "[N+;!R;X4](-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3])"), 
("[C10H21N]", "[N+;!R;X4](-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3])"),
("[C12H25N]", "[N+;!R;X4](-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3])"),
("[CH3P]", "[CH3]-[P+;!R;X4]"),    
("[C2H5P]", "[P+;!R;X4]-[CH2]-[CH3]"),    
("[C4H9P]", "[P+;!R;X4]-[CH2]-[CH2]-[CH2]-[CH3]"),     
("[C6H13P]", "[P+;!R;X4]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]"),    
("[C8H17P]", "[P+;!R;X4]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]"), 
("[C14H29P]", "[P+;!R;X4]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]"),  
("[CH3S]", "[S+;!R;X3](-[CH3])"),    
("[C2H5S]", "[S+;!R;X3](-[CH2]-[CH3])"),    
("[DABCO]", "[N+]-1-2-[C]-[C]-[N](-[C]-[C]-1)-[C]-[C]-2"),    
("[Quinc]", "[N+]-1-[C]-[C]-2-[C]-[C]-[C]-[C]-[C]-2-[C]-[C]-1"),
("[CL]", "[Cl-]"), 
("[BR]", "[Br-]"),    
("[BF4]", "[B-](-[F])(-[F])(-[F])-[F]"),    
("[PF6]", "[P-](-[F])(-[F])(-[F])(-[F])(-[F])-[F]"),    
("[AC]", "[C](=[O])(-[O-])-[C]"),     
("[14AC]", "[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C](-[O-])=[O]"),    
("[16AC]", "[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C](-[O-])=[O]"),     
("[TFA]", "[C](-[C](-[O-])=[O])(-[F])(-[F])-[F]"),    
("[NO3]", "[N+](-[O-])(=[O])-[O-]"),     
("[ClO4]", "[O-]-[Cl](~[O])(~[O])~[O]"),    
("[LA]", "[C]-[C](-[O])-[C](=[O])-[O-]"),    
("[SO3-PH1]", "c1c(ccc(c1)-[S](-[O-])(=[O])=[O])-[C]"),    
("[HSO4]", "[S](=[O])(-[O])(=[O])-[O-]"),    
("[SO3-O2O1]", "[C]-[O]-[C]-[C]-[O]-[S](-[O-])(=[O])=[O]"),     
("[SO3-O2O2]", "[C]-[C]-[O]-[C]-[C]-[O]-[S](-[O-])(=[O])=[O]"),    
("[SO3-O1]", "[S](=[O])(=[O])(-[O-])-[O]-[C]"),    
("[SO3-O2]", "[S](=[O])(=[O])(-[O-])-[O]-[C]-[C]"),     
("[SO3-O8]", "[S](=[O])(=[O])(-[O-])-[O]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]"),     
("[SO3-1]", "[S](=[O])(=[O])(-[O-])-[C]"),    
("[OTF]", "[S](=[O])(=[O])(-[C](-[F])(-[F])-[F])-[O-]"),     
("[C4F9SO3]", "[S](=[O])(=[O])(-[C](-[C](-[C](-[C](-[F])(-[F])-[F])(-[F])-[F])(-[F])-[F])(-[F])-[F])-[O-]"),    
("[SO3-O[2O]21]", "[S](=[O])(=[O])(-[O-])-[O]-[C]-[C]-[O]-[C]-[C]-[O]-[C]"),    
("[SCN]", "[S-]-[C]#[N]"),     
("[DCA]", "[N-](-[C]#[N])-[C]#[N]"),    
("[CCN3]", "[C-](-[C]#[N])(-[C]#[N])-[C]#[N]"),    
("[TCB]", "[B-](-[C]#[N])(-[C]#[N])(-[C]#[N])-[C]#[N]"),    
("[PO2-O1,O1]", "[P](-[O]-[C])(-[O]-[C])(=[O])-[O-]"),    
("[PO2-O2,O2]", "[O-]-[P](=[O])(-[O]-[C]-[C])-[O]-[C]-[C]"), 
("[PO2-O4,O4]", "[O-]-[P](=[O])(-[O]-[C]-[C]-[C]-[C])-[O]-[C]-[C]-[C]-[C]"),     
("[NTF2]", "[S](=[O])(=[O])(-[C](-[F])(-[F])-[F])-[N-]-[S](=[O])(=[O])-[C](-[F])(-[F])-[F]"),     
("[FAP]", "[P-](-[C](-[C](-[F])(-[F])-[F])(-[F])-[F])(-[C](-[C](-[F])(-[F])-[F])(-[F])-[F])(-[C](-[C](-[F])(-[F])-[F])(-[F])-[F])(-[F])(-[F])-[F]"),     
("[FSI]", "[N-](-[S](=[O])(=[O])-[F])-[S](=[O])(=[O])-[F]"),     
("[BOB]", "[B-]-1-2(-[O]-[C](-[C](-[O]-1)=[O])=[O])-[O]-[C](-[C](-[O]-2)=[O])=[O]"),     
("[TDI]", "[F]-[C](-[F])(-[F])-c1nc(c([n-]1)-[C]#[N])-[C]#[N]"),     
("[NPF2]", "[S](=[O])(=[O])(-[C](-[F])(-[C](-[F])(-[F])-[F])-[F])-[N-]-[S](=[O])(=[O])-[C](-[F])(-[C](-[F])(-[F])-[F])-[F]"),     
("[PO2-8I,8I]", "[C]-[C](-[C]-[C](-[C])(-[C])-[C])-[C]-[P](-[O-])(=[O])-[C]-[C](-[C])-[C]-[C](-[C])(-[C])-[C]"),     
("[PO2H-O1]", "[P](-[O-])(=[O])-[O]-[C]"),    
("[PO2H-O2]", "[P](-[O-])(=[O])-[O]-[C]-[C]"),
("[SBF6]", "F[Sb-](F)(F)(F)(F)F"),
("[TS]", "[O-]C(=O)c1ccccc1S"),
("[CPS]", "CCO[S](=O)(=O)C[C@]12CCC(CC1=O)C2(C)C")]

print(len(UNIFAC_SMARTS))
for i in UNIFAC_SMARTS:
    print(i)

# get the fragmentation scheme in the format necessary
fragmentation_scheme = {i+1: j[1] for i, j in enumerate(UNIFAC_SMARTS)}

# sort the fragmentation scheme according to the descriptors
pattern_descriptors = {}
for group_number, SMARTS in fragmentation_scheme.items():
    if type(SMARTS) is list:
        SMARTS = SMARTS[0]
    
    if SMARTS != "":
        pattern = fragmenter.get_mol_with_properties_from_SMARTS(SMARTS)
        
        pattern_descriptors[group_number] = [pattern.GetUnsignedProp('n_available_bonds') == 0,                                   (pattern.GetBoolProp('is_simple_atom_on_c') or pattern.GetBoolProp('is_simple_atom')),                                   pattern.GetUnsignedProp('n_atoms_defining_SMARTS'),  
                                  pattern.GetUnsignedProp('n_available_bonds') == 1, \
                                  fragmenter.get_heavy_atom_count(pattern) - pattern.GetUnsignedProp('n_carbons'), \
                                  pattern.GetBoolProp('has_atoms_in_ring'), \
                                  pattern.GetUnsignedProp('n_triple_bonds'), \
                                  pattern.GetUnsignedProp('n_double_bonds')]

sorted_pattern_descriptors = sorted(pattern_descriptors.items(), key=operator.itemgetter(1), reverse=True)
sorted_group_numbers = [i[0] for i in sorted_pattern_descriptors]


# In[ ]:


# second step: try to fragent all from the component    
structures_DB = []
with open('structures_DB.csv') as f:
    for line in f.readlines():
        structures_DB.append(CSV_to_info(line))
        
combined_fragmenter = fragmenter(fragmentation_scheme, 'combined', 20, function_to_choose_fragmentation, 1)
combined_fragmenter.fragmentation_scheme_order = sorted_group_numbers 
combined_fragmenter.n_max_fragmentations_to_find = 1    

combined_fragmenter_sorted_fragmented = []
right_size_for_combined_fragmenter  = []
print('####################################################################')
print('Fragmenting the structures database with the patterns sorted (combined algorithm)')
i_structure = 0
f_combined = open('structures_DB_combined_fragmentation_with_pattern_sorting_results.log','w+')
for inchikey, SMILES, pubchem_id, empty_fragmentation in structures_DB:  
 
    i_structure = i_structure + 1
    if i_structure % 4000 == 0:
        print('{:2.1f} .'.format((100.0 * i_structure) / len(structures_DB)), end=" ")
    
    fragmentation, success = combined_fragmenter.fragment(SMILES)
    if success:
        combined_fragmenter_sorted_fragmented.append(inchikey)
        
    n_heavy_atoms  = 0
    for sub_SMILES in SMILES.split("."):
        n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
    
    if n_heavy_atoms <= 20:
        right_size_for_combined_fragmenter.append(inchikey)
        log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {})
    else:
        if success:
            log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {})
        else:
            log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {}, 'Structure was skipped because it is larger than 20 atoms.')
          
  
f_combined.close()
    
print('')
print('N_structures(simple): ' + str(len(structures_DB)))
print('N_fragmented(simple): ' + str(len(combined_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(combined_fragmenter_sorted_fragmented)) / len(structures_DB)) + ")")
print('')
print('####################################################################')


# In[ ]:





# In[ ]:




