1. The chemical structure fragmentation of green solvents and solutes was under the environment that could be imported by the **my-rdkit-env.yml** file. Please refer to the **Environment** folder.
2. The scripts and examples for green solvents chemical structure fragmentation is under the folder **ILs and DESs**.
3. The scripts and examples for molecular solutes chemical structure fragmentation is under the folder **Solutes**.
4. As follows are the detailed instructions on how to do the chemical structure fragmentation taking the molecular solutes chemical structure fragmentation as an example:

    4.1 Under the **Solutes** folder, the **structures_DB.csv** file contains the SMILES for molecular solutes. 

    4.2 The **Fragmenter.ipynb** file or the **Fragmenter.py** file was the fragmentation script to conduct the fragmentation for molecular solutes representing by SMILES in the **structures_DB.csv** file. The preliminary fragmentation results were stored in the **structures_DB_combined_fragmentation_with_pattern_sorting_results.log** file.

    4.3 The **Manage_Solute_Fragmentation_Results.ipynb** file or the **Manage_Solute_Fragmentation_Results.py** file was the script to assemble the preliminary fragmentation results into fragmentation results of inputable format for the model. The fragmentation results of inputable format for the model was store in the **Solute_Fragmentation_Results.txt** file.

    4.4 The fragmentation results of inputable format for the model in the **Solute_Fragmentation_Results.txt** file was arranged in the **Solute_Fragmentation_Results.xlsx** file, and the fragmentation results in the sheet3 was the final fragmentation results of the solutes listed in the **structures_DB.csv** file.
