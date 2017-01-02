from Bio.KEGG import REST

#Get all pathways for organism (human in this case)
organism_pathways = REST.kegg_list("pathway","hsa").read()

#Get all reference pathways
reference_pathways = REST.kegg_list("pathway").read()

# For each organism pathway...
for org_line in organism_pathways.rstrip().split("\n"):

    #Get the pathway name and description
    org_entry, org_description = org_line.split("\t")
   
    #For each reference pathway...
    for ref_line in reference_pathways.rstrip().split("\n"):
        
        #Get the pathway name and description
        ref_entry, ref_description = ref_line.split("\t")

        #Check the current organism pathway description against the current reference pathway description, if it's the same...
        if org_description.find(ref_description) != -1:

           #Pull all EC #'s for reference pathway from KEGG using the reference pathway entry ('mapXXXX')
           org_ref_pathways_EC_nums = REST.kegg_link("ec",ref_entry).read()

           #Print pathway names and EC #'s
           print("All enzymes in human pathway "+org_entry+","+org_description+", which is reference pathway "+ref_entry+","+ref_description)
           print(org_ref_pathways_EC_nums)
        
