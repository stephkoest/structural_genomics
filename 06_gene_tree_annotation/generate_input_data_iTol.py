#!/usr/bin/env python3

from ete3 import Tree
import argparse
import pandas as pd
import os
import re
import collections

# This script generates iTol datasets to annotate a gene phylogeny with domain, phylum, class and taxonomy for prokaryotes (GTDB), the relevant taxonomy for eukaryotes and the taxonomy for Asgards
# It requires a newick tree file and metadata tables for archaea, bacteria, eukaryotes and Asgard archaea (Alpaca), as well as a list of leaf names for rooting the tree
# It uses Python library ETE3 to manipulate the tree and pandas to read the metadata tables
# It outputs a newick tree ('reformatted') and three iTol datasets: one for branch colours according to the species domain (Eukaryota, Archaea, Bacteria), one for new names and one for paralogs
# The script served as a basis to visualize single gene trees in iTol for the manuscript "Structure-based inference of eukaryotic complexity in Asgard archaea." by Koestlbacher et al., 2025

# Use the strings containing the taxonomies to figure out their overlap, which corresponds to the lowest taxonomic level shared by all
def find_largest_common_prefix(strings):
    "Across a list of strings, find the largest common string (present in all list items), from the left to the right"
    if not strings:
        return ""

    # Sort the list to bring the strings with the common prefix together
    strings.sort()

    # Consider the first and last strings after sorting
    first_str = strings[0]
    last_str = strings[-1]

    common_prefix = ""
    
    # Iterate through the characters of the first string
    for i in range(len(first_str)):
        # Check if the character is the same in the last string
        if first_str[i] == last_str[i]:
            common_prefix += first_str[i]
        else:
            break  # Stop if a difference is found

    return common_prefix

# Append a trailing slash to the directory name if it is not present
def checktrailingslash(directory_name):
    if directory_name.endswith("/") == False:
        directory_name += "/"
    return directory_name

# Read a file and return a list of lines without the newline character
def make_list_from_lines(input_file):
    with open(f"{input_file}", 'r') as infile:
        flines = infile.readlines()
        list_from_lines = [l.rstrip("\n") for l in flines]
    return list_from_lines

# Define the color scheme for each leaf prefix
color_scheme = collections.OrderedDict(
    {
    "Eukaryota": "#CC79A7",
    "Asgardarchaeota": "#009E73",
    "Archaea": "#0072B2",
    "Bacteria": "#D55E00"
    }
)

# Define Asgard nomenclature
asgard_shorttax_nomenclature = {
    "ALCG": "Atabeyarchaeia",
    "Atab": "Atabeyarchaeia",
    "Bald": "Baldrarchaeia",
    "Gefi": "Asgardarchaeia",
    "HeGe": "Gerdarchaeales",
    "HeHe": "Heimdallarchaeaceae",
    "HeHo": "Hodarchaeales",
    "HeKa": "Kariarchaeaceae",
    "Hela": "Helarchaeales",
    "HeNj": "Njordarchaeales",
    "Herm": "Hermodarchaeia",
    "Jord": "Jordarchaeia",
    "Loki": "Lokiarchaeales",
    "Odin": "Odinarchaeia",
    "Rana": "Ranarchaeia",
    "Sif": "Sifarchaeia",
    "Sifa": "Sifarchaeia",
    "Thor": "Thorarchaeia",
    "Wuko": "Wukongarchaeia"
}

# Shorten the leaf names to the first part of the name (before the first "/")
def simplify_leaf_names(tree):
    for leaf in tree:
        leaf.name = leaf.name.split("/")[0]

# Reroot the tree based on the leaves provided
def reroot_tree(tree, root_leaves):
    try: 
        if len(root_leaves) == 1:
            outgroup_node = tree.search_nodes(name=root_leaves[0])[0]
            tree.set_outgroup(outgroup_node)
        elif len(root_leaves) == 2:
            ancestor = tree.get_common_ancestor(root_leaves)
            tree.set_outgroup(ancestor)
    except AttributeError:
        midpoint = tree.get_midpoint_outgroup() 
        tree.set_outgroup(midpoint)

# Clean the tree from AlphaFold predicted structures that were added later and used as guidance (they were already represented)
def clean_tree(tree, leaves_to_remove=[]):
    """Remove branches corresponding to an AlphaFold predicted structure from the phylogeny"""
    if any(l.name.startswith("AF-") for l in tree.get_leaves()) or any(".pdb" in l.name for l in tree.get_leaves()) or len(leaves_to_remove) > 0:
        leaves_to_preserve = [l for l in tree.get_leaves() if l.name.startswith("AF-") == False and ".pdb" not in l.name and l not in leaves_to_remove]
        tree.prune(leaves_to_preserve, preserve_branch_length=True)

# Get the furthest descendants of a node
def get_furthest_descendants(node):
    if node.is_leaf() == False:
        child1, child2 = node.get_children()[0], node.get_children()[1]
        # get a random leaf of each child
        leaf1 = child1.get_leaf_names()[0]
        leaf2 = child2.get_leaf_names()[-1]
        return [leaf1, leaf2]
    else:
        return [node.name]

# Label the leaves with the domain, phylum, class and taxonomy for prokaryotes (GTDB), the relevant taxonomy for eukaryotes and the taxonomy for Asgards
def label_leaves(tree, arch_metadata, bact_metadata, euk_metadata, asg_metadata):
    """Assign labels to leaves (domain, phylum, class and taxonomy for prokaryotes (GTDB), the relevant taxonomy for eukaryotes), the taxonomy for Asgards"""
    unannotated_seqs = []
    for l in tree.get_leaves():
        l.domain = ""
        if l.name.startswith("Arch"):
            # Annotate archaea
            l.domain = "Archaea"
            l.phylum, l.clas = l.name.split("_")[1], l.name.split("_")[2]
            l.accession = re.search(r'_((GB|RS)_[^.]+\.\d+)_', l.name).group(1)
            if l.accession not in arch_metadata.index:
                print(f"Error: {l.accession} not found in metadata Archaea")
            else: 
                l.taxonomy = arch_metadata.loc[l.accession, "gtdb_taxonomy"]
                l.species = l.taxonomy.split("__")[-1]

        elif l.name.startswith("Bact"):
            # Annotate bacteria
            l.domain = "Bacteria"
            l.phylum, l.clas = l.name.split("_")[1], l.name.split("_")[2]
            l.accession = re.search(r'_((GB|RS)_[^.]+\.\d+)_', l.name).group(1)
            if l.accession not in bact_metadata.index:
                print(f"Error: {l.accession} not found in metadata Bacteria")
            else:
                l.taxonomy = bact_metadata.loc[l.accession, "gtdb_taxonomy"]
                l.species = l.taxonomy.split("__")[-1]
        
        else:
            # Annotate eukaryotes and Asgards
            if re.search(r'^(\D{6})(\d{6})', l.name):
                l_species_acronym = re.search(r'^(\D{6})(\d{6})', l.name).group(1)
                if l_species_acronym in euk_metadata.index:
                    l.species = euk_metadata.loc[l_species_acronym, "Scientific name"]
                    l.domain = "Eukaryota"
                    l.clade = euk_metadata.loc[l_species_acronym, "relevant taxonomy"]
                else:
                    print(f"{l.name} has eukaryote-like identifier but is not found in eukaryotic metadata dataframe")
            elif re.search(r'^(\S+)_(\d+)$', l.name):
                l_species_acronym = re.search(r'^(\S+)_(\d+)$', l.name).group(1)
                if l_species_acronym in asg_metadata.index:
                    l.accession = l_species_acronym
                    shorttaxname = asg_metadata.loc[l_species_acronym, "ShortTaxName"]
                    l.species = f'{shorttaxname} {asg_metadata.loc[l_species_acronym, "Strain"]}'
                    l.domain = "Asgardarchaeota"
                    l.clade = asgard_shorttax_nomenclature[shorttaxname] if shorttaxname in asgard_shorttax_nomenclature else ""
                else:         
                    print(f"{l.name} has an unmappable identifier ({l_species_acronym})")
        # List and remove leaves that have not been assigned a domain yet
        if l.domain == "":
            unannotated_seqs.append(l)
    clean_tree(tree, unannotated_seqs)


# Label the internal nodes with the protein name and the taxonomy for prokaryotes (GTDB)
def label_internal_nodes(tree):
    """Assign labels to internal nodes: protein names corresponding to the (orthologous) groups and the clade/taxon name and rank for prokaryotic gtdb clades only"""
    hierarchy = collections.OrderedDict({'s':'species', 'g':'genus', 'f':'family', 'o':'order', 'c':'class', 'p':'phylum', 'd':'domain'})
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            n.domain = ""
            n.clade = ""
            n.rank = ""
            leaves = n.get_leaves()
            # Then get the domain, an the clade and rank according to the leaf taxonomies (prokaryotic only)
            leaves_domain = [l.domain for l in leaves]
            if all(item == leaves_domain[0] for item in leaves_domain):
                n.domain = leaves_domain[0]
                # Check if in GTDB based on domain assignment; if so find the lowest taxonomic level and its rank
                if n.domain == "Archaea" or n.domain == "Bacteria":
                    leaves_taxonomies = [l.taxonomy for l in leaves]
                    n_taxonomy_full = find_largest_common_prefix(leaves_taxonomies)
                    ## If not shared up to the species level, the pattern will end with the rank and "__", which should be removed
                    if n_taxonomy_full.endswith("__"):
                        n_taxonomy_full = n_taxonomy_full[0:-4]
                    # Select the lowest level
                    n_taxonomy_specific = n_taxonomy_full.split(";")[-1]
                    n.rank, n.clade = n_taxonomy_specific.split("__")
                    n.rank = hierarchy[n.rank]
                # If not in GTDB (Eukaryota, Asgardarchaeota), check if the "clade" is the same for all leaves - if so annotate as such
                elif all(l.clade == leaves[0].clade for l in leaves): 
                    n.clade = leaves[0].clade
                    n.rank = "NA"         
    # Visit the nodes once more to label unannotated nodes based on parent and tips
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            if n.domain == "":
                try: 
                    p_domain = n.up.domain
                    l_domains = [c.domain for c in n.get_leaves()]
                    # Check if there's at least one match between parent and children (the other children might have another domain annotation)
                    if any(l_domain == p_domain for l_domain in l_domains):
                        n.domain = p_domain
                except AttributeError:
                    pass
                
# Generate the iTol dataset for branch colours according to the species domain (Eukaryota, Archaea, Bacteria)
def generate_itol_dataset_branch_colours(tree):
    """Generate iTol dataset from a tree and a root prefix"""
    # create the iTol dataset in comma-separated format
    dataset = []
    # Find monophyletic groups of bacteria, archaea and eukaryotes and colour them (the ancestor and the entire clade emanating from it)
    for dom in color_scheme.keys():
        for n in tree.get_monophyletic(values=[dom], target_attr="domain"):
            descendants = get_furthest_descendants(n)
            dataset.append(f"{'|'.join(descendants)},branch,clade,{color_scheme[dom]},1,normal")
    # Also colour those internal nodes that aren't monophyletic with regard to the domain, but that do have an annotated domain (test)
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            for dom in color_scheme.keys():
                if n.domain == dom:
                    if n.check_monophyly(values=[dom], target_attr="domain") == False:
                        descendants = get_furthest_descendants(n)
                        dataset.append(f"{'|'.join(descendants)},branch,clade,{color_scheme[dom]},1,normal")
    return dataset

# Generate the iTol dataset for new names
def generate_itol_dataset_new_names(tree):
    """Generate iTol dataset that renames the taxa"""
    dataset = []
    for l in tree.get_leaves():
        if l.name.startswith(("Bact", "Arch")):
            dataset.append(f"{l.name},{l.phylum}_{l.species.replace(' ', '_')}")
        else: 
            dataset.append(f"{l.name},{l.clade}_{l.species.replace(' ', '_')}")
    # Label the internal nodes with the clade - if available
    for n in tree.iter_descendants("preorder"):
        if n.is_leaf() == False:
            descendants = get_furthest_descendants(n)
            if n.clade != "":
                dataset.append(f"{'|'.join(descendants)},{n.clade}")
    return dataset

# Generate a list of colours for the paralog leaves
def get_color_list(num_colors):
    import colorsys
    # Generate equally spaced hues
    hues = [i / num_colors for i in range(num_colors)]
    # Convert hues to RGB values
    rgb_values = [colorsys.hsv_to_rgb(hue, 1, 1) for hue in hues]
    # Convert RGB values to hex codes
    hex_list = ['#%02x%02x%02x' % tuple(int(round(255*val)) for val in color) for color in rgb_values]
    return hex_list

# Generate the iTol dataset for paralogs
def generate_itol_dataset_paralogs(tree, arch_metadata, bact_metadata, asg_metadata):
    """Add coloured symbols to paralog leaves"""
    dataset = []
    paralog_taxa = []
    # First get the taxa that have paralogs in the tree
    for l in tree.get_leaves():
        if l.domain == "Archaea" and l.accession in arch_metadata.index:
            if len(tree.search_nodes(domain="Archaea", accession = l.accession)) > 1 and l.name not in paralog_taxa:
                paralog_taxa.append(l.accession)
        elif l.domain == "Bacteria" and l.accession in bact_metadata.index:
            if len(tree.search_nodes(domain="Bacteria", accession = l.accession)) > 1 and l.name not in paralog_taxa:
                paralog_taxa.append(l.accession)
        elif (l.domain == "Asgardarchaeota") and l.accession in asg_metadata.index:
            if len(tree.search_nodes(domain="Asgardarchaeota", accession = l.accession)) > 1 and l.name not in paralog_taxa:
                paralog_taxa.append(l.accession)
    paralog_taxa = list(set(paralog_taxa))
    paralog_taxa.sort()
    # Add star symbols for each taxon with duplicates
    hex_range = get_color_list(len(paralog_taxa))
    for i, taxon in enumerate(paralog_taxa):
        taxon_nodes = tree.search_nodes(accession = taxon)
        for paralog in taxon_nodes:
            dataset.append(f"{paralog.name},3,0.5,{hex_range[i]},1,0.8")
    return(dataset)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate iTol datasets to annotate a gene phylogeny')
    parser.add_argument('-t', metavar='tree_path', type=str, help='Path to the input newick tree file')
    parser.add_argument('-ma', metavar='archaea', type=str, help='Path to metadata table of archaeal lineages', default="phylogeny_annotation_files_ITOL/ar53_metadata_r207.qscore.family_representative.csv")
    parser.add_argument('-mb', metavar='bacteria', type=str, help='Path to metadata table of bacterial lineages', default="phylogeny_annotation_files_ITOL/bac120_metadata_r207.qscore.family_representative.csv")
    parser.add_argument('-me', metavar='eukaryota', type=str, help='Path to metadata table for eukaryotes', default="phylogeny_annotation_files_ITOL/Euk5FinalSet.adjust.busco.euk5_tree_abbrev.csv")
    parser.add_argument('-ms', metavar='asgards', type=str, help='Path to metadata table for Asgard archaea (Alpaca)', default="phylogeny_annotation_files_ITOL/Asgard_DB_230420.tsv")
    parser.add_argument('-o', metavar='output_dir', type=str, help='output directory - current working directory if not specified')
    parser.add_argument('-r', metavar='root_leaves', nargs='+', type=str, help='List of leaf names for rooting the tree')
    args = parser.parse_args()

    # Get the output directory and make it if it does not exist yet
    if args.o == None:
        outdir = checktrailingslash(os.getcwd())
    else:
        outdir = checktrailingslash(args.o)
        try:
            os.mkdir(outdir)
        except FileExistsError:
            pass

    # Get the basename of the input tree
    tree_basename = os.path.basename(args.t)
    
    # Load the tree from file
    tree = Tree(args.t)

    # Simplify leaf names
    simplify_leaf_names(tree)

    # Reroot the tree
    reroot_tree(tree, args.r)
    
    # Load species metadata
    archaea_metadata = pd.read_csv(args.ma, index_col="accession")
    bacteria_metadata = pd.read_csv(args.mb, index_col="accession")
    eukaryota_metadata = pd.read_csv(args.me, index_col="Abbreviation")
    asgard_metadata = pd.read_csv(args.ms, index_col="Locus tag (Strain no underscore)", sep="\t")
    # Remove duplicates in Asgard metadata, keep the first entry
    asgard_metadata = asgard_metadata[~asgard_metadata.index.duplicated(keep='first')]
    
    # Remove sequences from AlphaFold structures from the tree
    clean_tree(tree)
    
    # Label the leaves according to their domain and lower taxonomy or proteins + remove non-labelled sequences with clean_tree
    label_leaves(tree, archaea_metadata, bacteria_metadata, eukaryota_metadata, asgard_metadata)

    # Write new tree to newick
    tree.write(outfile=f"{outdir}{tree_basename}.reformatted", format=0)

    # Label the internal nodes: protein name and taxonomy (the latter for prokaryotes only)
    label_internal_nodes(tree)

    # Generate and write the iTol dataset to branch colours according to the species domain (Eukaryota, Archaea, Bacteria)
    dataset_branch_colours = generate_itol_dataset_branch_colours(tree)
    with open(f"{outdir}{tree_basename}.reformatted.iTOL_domain.dataset.txt", "w") as file:
        file.write("DATASET_STYLE\n")
        file.write("SEPARATOR COMMA\n")
        file.write("DATASET_LABEL,Domain\n")
        file.write("COLOR,#ffff00\n")
        file.write("LEGEND_TITLE,Domain\n")
        file.write("LEGEND_POSITION_X,100\n")
        file.write("LEGEND_POSITION_Y,100\n")
        file.write("LEGEND_HORIZONTAL,0\n")
        file.write(f"LEGEND_SHAPES{',1' * len(color_scheme)}\n")
        file.write(f"LEGEND_COLORS,{','.join(color_scheme.values())}\n")
        file.write(f"LEGEND_LABELS,{','.join(color_scheme.keys())}\n")
        file.write(f"LEGEND_SHAPE_SCALES{',1' * len(color_scheme)}\n")
        file.write("DATA\n")
        for item in dataset_branch_colours:
            file.write(f"{item}\n")    

    # Generate and write the iTol dataset for new names
    dataset_labels = generate_itol_dataset_new_names(tree)
    with open(f"{outdir}{tree_basename}.reformatted.iTOL_labels.dataset.txt", "w") as file:
        file.write("LABELS\n")
        file.write("SEPARATOR COMMA\n")
        file.write("DATA\n")
        for item in dataset_labels:
            file.write(f"{item}\n")
    
    # Generate and write the iTol dataset for paralogs 
    dataset_paralogs = generate_itol_dataset_paralogs(tree, archaea_metadata, bacteria_metadata, asgard_metadata)
    with open(f"{outdir}{tree_basename}.reformatted.iTOL_paralogshapes.dataset.txt", "w") as file:
        file.write("DATASET_SYMBOL\n")
        file.write("SEPARATOR COMMA\n")
        file.write("DATASET_LABEL,Paralogs prokaryotes\n")
        file.write("COLOR,#AC3A6D\n")
        file.write("LEGEND_TITLE,Paralogs prokaryotes\nLEGEND_POSITION_X,80\nLEGEND_POSITION_Y,80\nLEGEND_HORIZONTAL,0\nLEGEND_SHAPES,3\nLEGEND_COLORS,#AC3A6D\nLEGEND_LABELS,paralog\nLEGEND_SHAPE_SCALES,1\nLEGEND_SHAPE_INVERT,0\n")
        file.write("MAXIMUM_SIZE,10\n")
        file.write("#GRADIENT_FILL,1\n")
        file.write("DATA\n")
        for item in dataset_paralogs:
            file.write(f"{item}\n")
    

                   













