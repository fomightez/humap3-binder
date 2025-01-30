# Make a list of all accession identifiers shared in both version 2.0 and 
# version 3.0 of data. And then use those to go through and collect identifiers 
# for 'companion' proteins that occur in all of the associated version 2.0
# complexes, but then go disappearing entirely from the related set of version 
# 3.0 complexes.
# While at it collect, the identifiers used to collect each member of 'the 
# ubiquitous disappearing set'. These will be useful candidate proxies to
# investigate where the disappearing proteins would be expected in version 3.0
# of the data. Plus, while at it also collect the hu.MAP2.0 complex identifiers
# involved.
# To be strict in order to possibly find reasonable candidates where 
# 'actual-associated' proteins seen in hu.MAP 2.0 data complexes go missing in 
# related hu.MAP 3.0 data complexes, I stipulate the number of complexes needed 
# to be more than one to qualify as 'ubiquitous' in hu.MAP 2.0 complexes.
# In other words, to be 'ubiquitous', has to disappear from a multiple 
# related complexes. These are called `ubiquitous_disappearing_list/ubiquitous_disappearing_set/disappearing_identified_with_dict/disappeared_total_complexes_dict`. 
# (Note that turns out to highlight about 10% [or less if you subtract the 1.5% 
# that totally disappear in version 3] where things shift around slightly 
# between the two versions is what is specifically being identified by all but 
# a small percent of that group.)
# Note this ended up being 2055 identifiers, about one-fifth of the identifiers 
# in the hu.MAP 2.0 data, and so I added another level where I only collected
# those that weren't in any complexes at all among the hu.MAP 3.0 data. These
# I called `gone_entirely_in_3_list/gone_entirely_in_3_set/gone_entirely_in_3_identified_via_dict/gone_entirely_in_3_total_complexes_dict`. Can think of this level removing those 
#  aren't in particular 'related' complexes but still seen in v3 at least.

#First, set up getting the data in memory
import os
import pandas as pd
csv_file_raw_data2 = "humap2_complexes_20200809InOrderMatched.csv"
if not os.path.isfile(csv_file_raw_data2):
    os.system("curl -OL -s https://raw.githubusercontent.com/"\
        "fomightez/humap2-binder/refs/heads/main/additional_nbs/"\
        f"standardizing_initial_data/{csv_file_raw_data2}")
csv_file_raw_data = "hu.MAP3.0_complexes_wConfidenceScores_total15326_wGenenames_20240922InOrderMatched.csv"
if not os.path.isfile(csv_file_raw_data):
    os.system("curl -OL -s https://raw.githubusercontent.com/"\
        "fomightez/humap3-binder/refs/heads/main/additional_nbs/"\
        f"standardizing_initial_data/{csv_file_raw_data}")
CSV2df_script_needed = "complexes_rawCSV_to_df.py"
if not os.path.isfile(CSV2df_script_needed):
    os.system("curl -OL -s https://raw.githubusercontent.com/"\
        "fomightez/structurework/refs/heads/master/"\
        f"humap3-utilities/{CSV2df_script_needed}")
lud_script_needed = "make_lookup_table_for_extra_info4complexes.py"
if not os.path.isfile(lud_script_needed):
    os.system("curl -OL -s https://raw.githubusercontent.com/"\
        "fomightez/structurework/refs/heads/master/"\
        f"humap3-utilities/{lud_script_needed}")
os.system(f"uv run complexes_rawCSV_to_df.py {csv_file_raw_data2}")
rd_df2 = pd.read_pickle('raw_complexes_pickled_df.pkl')
print("\n")
os.system(f"uv run complexes_rawCSV_to_df.py {csv_file_raw_data}")
print("\n") # add exta line break to set up for next `print()`; otherwise prints
# far over on right
rd_df3 = pd.read_pickle('raw_complexes_pickled_df.pkl')

# Now with data in memory, to get all the accessions shared by both first I need 
# to split out all and get intersection.
d2_4_expanding = rd_df2.copy()
d3_4_expanding = rd_df3.copy()
d2_4_expanding['Uniprot_ACCs'] = d2_4_expanding['Uniprot_ACCs'].str.split()
d3_4_expanding['Uniprot_ACCs'] = d3_4_expanding['Uniprot_ACCs'].str.split()
v2data_info_as_df = d2_4_expanding.explode(['Uniprot_ACCs']).copy()
v3data_info_as_df = d3_4_expanding.explode(['Uniprot_ACCs']).copy()
v2data_info_as_df = v2data_info_as_df.reset_index(drop=True)
v3data_info_as_df = v3data_info_as_df.reset_index(drop=True)
accessions_shared_by_both_2and3 = set(v2data_info_as_df['Uniprot_ACCs']).intersection(set(v3data_info_as_df['Uniprot_ACCs']))
disappearing_identified_with_dict = {} # start the dictionary storing the 
# accession identifiers being examined when disappearing ones noted; keys will 
# be the disappearing ones.
# Actually, the keys for `disappearing_identified_with_dict` will be the same
# identifiers(accs) as I was planning to collect as list and then make into a 
# set to drop the duplicates as 'the ubiquitous disappearing set' , so I can just 
# use that to assign that later.
v2_complexes_with_disappearing_dict = {} # start the dictionary storing the 
#hu.MAP2.0 complex ids related to the disappearing ones; keys will be the disappearing ones
disappeared_total_complexes_dict = {} # start the dictionary storing the total 
#complexes in hu.MAP 2.0 for each member of the 
# 'the ubiquitous disappearing set'; the keys will also be the disappearing ones

# The iterating takes a fair amount of time so I went ahead and saved the pickled verisons of the dictionaries that should be made.
# Detect the main pickled dict file being present as a way to tell if really want to do the iterating or just read the pickled files
main_pickled_dict_file_to_check_for = 'disappearing_identified_with_dict.pkl'
import pickle 
def unpickle_dict(file_name):
    with open(file_name, "rb") as f:
        return pickle.load(f)
if os.path.exists(main_pickled_dict_file_to_check_for):
    disappearing_identified_with_dict = unpickle_dict(main_pickled_dict_file_to_check_for)
    v2_complexes_with_disappearing_dict = unpickle_dict("v2_complexes_with_disappearing_dict.pkl")
    disappeared_total_complexes_dict = unpickle_dict("disappeared_total_complexes_dict.pkl")
    gone_entirely_in_3_identified_via_dict = unpickle_dict("gone_entirely_in_3_identified_via_dict.pkl")
    v2_complexes_dict_for_gone_in_3 = unpickle_dict("v2_complexes_dict_for_gone_in_3.pkl")
    gone_entirely_in_3_total_complexes_dict = unpickle_dict("gone_entirely_in_3_total_complexes_dict.pkl")
else:
    for search_term in accessions_shared_by_both_2and3:
        pattern = fr'\b{search_term}\b' # Create a regex pattern with word boundaries
        # Set up and getting list of all associated in both version 2 and version 3 data and make set of missing from version 3 (based on `Highlight_differences_between_hu.MAP2_and_hu.MAP3_data.ipynb` and `Highlight_differences_between_hu.MAP2_and_hu.MAP3_data.ipynb` and sections 'Show all complexes that protein is in with extra information' and 'Show all proteins in related complexes with details added from Uniprot' in the notebook `Working_with_hu.MAP3_data_with_Python_in_Jupyter_Basics.ipynb`)
        d2_rows_with_term_df = rd_df2[rd_df2['Uniprot_ACCs'].str.contains(pattern, case=False, regex=True) | rd_df2['genenames'].str.contains(pattern, case=False, regex=True)].copy()
        d3_rows_with_term_df = rd_df3[rd_df3['Uniprot_ACCs'].str.contains(pattern, case=False, regex=True) | rd_df3['genenames'].str.contains(pattern, case=False, regex=True)].copy()
        d2_rows_with_term_df['Uniprot_ACCs'] = d2_rows_with_term_df['Uniprot_ACCs'].str.split()
        d3_rows_with_term_df['Uniprot_ACCs'] = d3_rows_with_term_df['Uniprot_ACCs'].str.split()
        d2_rows_with_term_df['genenames'] = d2_rows_with_term_df['genenames'].str.split()
        d3_rows_with_term_df['genenames'] = d3_rows_with_term_df['genenames'].str.split()
        df2_expanded = d2_rows_with_term_df.explode(['Uniprot_ACCs', 'genenames']).copy()
        df3_expanded = d3_rows_with_term_df.explode(['Uniprot_ACCs', 'genenames']).copy()
        df2_expanded = df2_expanded.reset_index(drop=True)
        df3_expanded = df3_expanded.reset_index(drop=True)
        df2_and_3_accs = set(df2_expanded['Uniprot_ACCs'].to_list() + df3_expanded['Uniprot_ACCs'].to_list()) # needed here only to make lookup table , or mock version
        #%run -i make_lookup_table_for_extra_info4complexes.py # This adds a lot of time and the additional information is currently not used; uncomment this line if you later decide you need that added information to explore more
        if 'lookup_dict' not in locals():
            # make a mock dictionary lookup_dict to use since not making lookup table in
            # the interest of saving time since that information presently not used.
            # Just uncomment the line above that is 
            # `%run -i make_lookup_table_for_extra_info4complexes.py` in order to restore
            # and this section will get skipped. 
            lookup_dict = {}
            for acc in df2_and_3_accs:
                mock_data_string = "LOOKUP TABLE WITH THIS INFORMATION NOT MADE. SEE CODE AND UNCOMMENT A LINE TO GET ACTUAL DATA."
                lookup_dict[acc] = {'protein_name':mock_data_string, 'disease': mock_data_string, 'synonyms': mock_data_string}
        pn_dict = {k: v['protein_name'] for k, v in lookup_dict.items()}
        disease_dict = {k: v['disease'] for k, v in lookup_dict.items()}
        synonyms_dict = {k: v['synonyms'] for k, v in lookup_dict.items()}
        df3_expanded['synonyms'] = df3_expanded['Uniprot_ACCs'].map(synonyms_dict)
        df3_expanded['protein_name'] = df3_expanded['Uniprot_ACCs'].map(pn_dict)
        df3_expanded['disease'] = df3_expanded['Uniprot_ACCs'].map(disease_dict)
        conf_val2text_dict = {
            1: 'Extremely High',
            2: 'Very High',
            3: 'High',
            4: 'Moderate High',
            5: 'Medium High',
            6: 'Medium'
        }
        df3_expanded['ComplexConfidence'] = df3_expanded['ComplexConfidence'].map(conf_val2text_dict)
        base_uniprot_url = 'https://www.uniprot.org/uniprotkb/'
        df3_expanded = df3_expanded.assign(Link=base_uniprot_url + df3_expanded['Uniprot_ACCs'])
        df2_expanded['synonyms'] = df2_expanded['Uniprot_ACCs'].map(synonyms_dict)
        df2_expanded['protein_name'] = df2_expanded['Uniprot_ACCs'].map(pn_dict)
        df2_expanded['disease'] = df2_expanded['Uniprot_ACCs'].map(disease_dict)
        conf_val2text_dict = {
            1: 'Extremely High',
            2: 'Very High',
            3: 'High',
            4: 'Moderate High',
            5: 'Medium High',
            6: 'Medium'
        }
        df2_expanded['Confidence'] = df2_expanded['Confidence'].map(conf_val2text_dict)
        base_uniprot_url = 'https://www.uniprot.org/uniprotkb/'
        df2_expanded = df2_expanded.assign(Link=base_uniprot_url + df2_expanded['Uniprot_ACCs'])
        accs_in2 = df2_expanded['Uniprot_ACCs'].to_list()
        accs_in3 = df3_expanded['Uniprot_ACCs'].to_list()
        missing_from_datav3 = list(set(accs_in2) - set(accs_in3))

        # Next make a list of those associated that occur in all complexes from data version 2 (based on `find_majority_uniprot_accs()` function in `Highlight_differences_between_hu.MAP2_and_hu.MAP3_data.ipynb`, but increase to not just be majority, but if occur in all associated complexes)
        # To be classified as 'ubiquitous' here, the number of associated 
        # complexes has to be more than one. May be too strict but was seeing 
        # there was a lot without that constraint, and so probably reasonable
        # for now.
        def find_ubiquitous_uniprot_accs(df):
            # Group by HuMAP2_ID and get unique Uniprot_ACCs for each complex
            accs_per_complex_id_series = df.groupby('HuMAP2_ID')['Uniprot_ACCs'].unique()
            
            # Count total number of unique complexes
            total_complexes = len(accs_per_complex_id_series)
            
            # Flatten and count occurrences of each Uniprot ACC across complexes
            all_accs = df['Uniprot_ACCs'].tolist()
            acc_counts = pd.Series(all_accs).value_counts()
            
            # Find ACCs that appear in ALL of the complexes
            if total_complexes > 1: # To be 'ubiquitous' the number of associated complexes has to be more than one. May be too strict but there was a lot without it, and so probably not at this time.
                ubiquitous_accs = acc_counts[acc_counts >= 1.00 * total_complexes].index.tolist()
            else:
                ubiquitous_accs = []
            return ubiquitous_accs
        if missing_from_datav3: # no point in continuing on if nothing missing
            ubiquitous_identifiers = find_ubiquitous_uniprot_accs(df2_expanded)

        # Function to check shared items
        def get_shared_items(set1, set2):
            return set1.intersection(set2)
        # Those that disappear from ubiquitous-occuring group will be the 
        # intersection of the two lists `missing_from_datav3` and 
        # `ubiquitous_identifiers` - when identify ones, record what 
        # identifier was the `search_term` was used to find this. Only collect 
        # those that have no 'weird 'notes in the genename.
        if missing_from_datav3 and ubiquitous_identifiers:
            ubiquitous_disappeared_identifiers = get_shared_items(set(missing_from_datav3),set(ubiquitous_identifiers))
            for x in ubiquitous_disappeared_identifiers:
                if x in disappearing_identified_with_dict:
                    disappearing_identified_with_dict[x].append(search_term)
                    # the next conditional will add the hu.MAP2.0 complex identifiers if that information concerning those complexes not yet added.
                    for i in df2_expanded.groupby('HuMAP2_ID')['Uniprot_ACCs'].unique().index.to_list():
                        if i not in v2_complexes_with_disappearing_dict[x]:
                            v2_complexes_with_disappearing_dict[x].append(i)
                    ''' THIS BELOW WAS WHEN I WAS JUST ADDING COMPLEX SIZES AS COME ALONG BUT IF RELATED ONES BOTH HAD '2' it would only add one of the two values, so probably not clearest way
                    # the next conditional will add the total complex size into a list if not already present (added later to see if there are multiple total number of assiciated complexes)
                    if len(df2_expanded.groupby('HuMAP2_ID')['Uniprot_ACCs'].unique()) not in disappeared_total_complexes_dict[x]:
                        disappeared_total_complexes_dict[x].append(len(df2_expanded.groupby('HuMAP2_ID')['Uniprot_ACCs'].unique()))
                    '''
                else:
                    disappearing_identified_with_dict[x] = [search_term]
                    v2_complexes_with_disappearing_dict[x] = df2_expanded.groupby('HuMAP2_ID')['Uniprot_ACCs'].unique().index.to_list()

    disappeared_total_complexes_dict = {k:len(cl) for k,cl in v2_complexes_with_disappearing_dict.items()}
    # Just ran this and it takes 12min 36s typically to get to this point iterating on the 9574 UniProt accension identifiers that occur BOTH in version 2 and 3 of the data, `accessions_shared_by_both_2and3`, to find the 2055 identifiers,` i.e., len(disappearing_identified_with_dict.keys())`, (out of 10058 total in hu.MAP 2.0 data) for the associated proteins that occur in all the associated complexes in vesion 2.0, but disappear in version 3.0 of the related complexes data. (Note it was 4303 identifiers found before I stipulated the number of complexes needed to be more than one to qualify as 'ubiquitous' in hu.MAP 2.0 complexes.)




    # There is a level above this where selection for disappeared is more extreme because only consider those gone entirely from version 3 data, not just gone in what seem to be related complexes based on what was seen in hu.MAP 2.0 complex members (TO MAKE THIS MORE
    # STRICT LEVEL, I CAN USE the `disappearing_identified_with_dict`/'ubiquitous_disappearing..'-related stuff and filer out any that are in hu.MAP 3.0 data, which means they aren't in particular 'related' complexes but still seen at least; the first two dictionary comprehensions take a bit of time to run)
    gone_entirely_in_3_identified_via_dict = {k:v for k,v in disappearing_identified_with_dict.items() if k not in set(v3data_info_as_df['Uniprot_ACCs'])}
    v2_complexes_dict_for_gone_in_3 = {k:v for k,v in v2_complexes_with_disappearing_dict.items() if k not in set(v3data_info_as_df['Uniprot_ACCs'])}
    gone_entirely_in_3_total_complexes_dict = {k:v for k,v in disappeared_total_complexes_dict.items() if k not in set(v3data_info_as_df['Uniprot_ACCs'])}

    # Save all those dictionaries with pickling so can develop remainder after without needing to run all those again
    #Serialize the data; use of `with` based on https://stackoverflow.com/a/20101064/8508004
    def pickle_dict(d,file_name):
        with open(file_name, "wb") as f:
            pickle.dump(d, f)
    pickle_dict(disappearing_identified_with_dict, "disappearing_identified_with_dict.pkl")
    pickle_dict(v2_complexes_with_disappearing_dict, "v2_complexes_with_disappearing_dict.pkl")
    pickle_dict(disappeared_total_complexes_dict, "disappeared_total_complexes_dict.pkl")

    pickle_dict(gone_entirely_in_3_identified_via_dict, "gone_entirely_in_3_identified_via_dict.pkl")
    pickle_dict(v2_complexes_dict_for_gone_in_3, "v2_complexes_dict_for_gone_in_3.pkl")
    pickle_dict(gone_entirely_in_3_total_complexes_dict, "gone_entirely_in_3_total_complexes_dict.pkl")