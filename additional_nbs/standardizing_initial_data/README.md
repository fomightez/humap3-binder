# Fixing and balancing the author-provided data to make a standardized CSV

This is 'Techincal'-related topic and not for you if you are just looking to explore complexes...

Dealing with the author-provided CSV file, I had noticed it was 'untidy' in that it didn't adhere to one variable per column and row. Pandas `explode()` is able to handle that but it turns out the identifiers between the two columns,  `Uniprot_ACCs` & `genenames`, were not in-order matched as I had initially through based on cursory examination of some manageable examples.
Looking into that 'untidy' data further I noted that in addition to data multiple identifiers per row, not in consistent order between the two columns that there's also mixed data types with 'mixed delimiter' issues & some inconsistencies.
The Jupyter notebook, and associted files, here help address that and make an in-order matched and balanced CSV file.   
See [the notebook itself 'Jupyter Notebook Standardizing Identifier Order In Adjacent Columns in hu.MAP3-provided CSV'](Standardizing_identifier_order_in_humap3-provided_csv.ipynb) for details and processing to standardize the CSV contents better.  

Note, for the shareable Jupyter sessions to make it easy to work through the Jupyter notebooks in the series, I already was not able to get the data directly from the author's site because MyBinder.org-provided sessions do not seem to allow access. I think a port is blocked. And so I had already copied the current version of the author-provided data to a resource MyBinder could access. So now the concern will be keeping it consistent with what the authors provide in case they update it. I did though leave this notebook here in a way that it can easily be run if the source data is extended or improved by the authors.

Quick-Link -- Related Notebook:
- [(Technical) Jupyter Notebook Standardizing Identifier Order In Adjacent Columns in hu.MAP3-provided CSV](Standardizing_identifier_order_in_humap3-provided_csv.ipynb)