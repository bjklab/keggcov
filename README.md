# keggcov

The *keggcov* package leverages the KEGG metabolic pathway database <https://www.genome.jp/kegg/pathway.html> to estimate a variance-covariance structure for sets of KOs identified by high-throughput sequencing. The output matrix can be interpreted through its (1) Diagonal Elements (Variance), which represent the functional versatility of the KO (i.e., a high variance on the diagonal means the KO is involved in many different pathways (a "hub" gene); a low variance means the KO is specialized to only one or two specific pathways); and through its (2) Off-Diagonal Elements (Covariance), which prepresent pathway partnerships. Positive high covariance indicates that the KOs are "pathway partners." They almost always appear together (e.g., subunits of the same complex or sequential enzymes in a linear pathway). Zeros indicated that the KOs share no biological processes (e.g., a DNA replication gene vs. a glycolysis gene).  

The package includes includes:  
+ a function that takes as input a set of KOs and outputs a variance-covariance matrix 

The *keggcov* package has been used to create data and/or figures for the following publications: .

