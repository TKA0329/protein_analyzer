# Protein Properties Analyzer
Link to Streamlit Web App: https://proteinanalyzer-tfdxptasphzhddkxctm3ed.streamlit.app/

## How it works

Put your protein sequences in a CSV file in the following format, with `sequence` as the header:

```csv
sequence
MPYEKHVEQTVVEKTEQGGSGGSYRHQTEAEAEKIRRELEKQGGGGSGGGGS
RQGGGAPAGGNIGGGQPQGGWGQPQQPQGGNQFSGGAQSRPQ
MKTQRDGHSLGRWSLVLLLLGLVMPLAIIAQVLSYKEAVL
```

Upload it to the web app and it will return a table with the following properties for each sequence: molecular weight, isoelectric point (pI), GRAVY score, aromaticity, instability index, net charge at pH 7, and secondary structure fractions.

Useful for processing large numbers of sequences at once without having to submit them one by one to tools like ProtParam.

