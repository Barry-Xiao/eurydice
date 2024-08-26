feces = {
    2: {'small': 4,
        'small_groups': [
            'p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae',
            'p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae',
            'p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Family_XI',
            'p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae',
            'p__Firmicutes;c__Other;o__Other;f__Other',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae',
            'p__Bacteroidota;c__Other;o__Other;f__Other',
            'p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae',
            'p__Actinobacteriota;c__Other;o__Other;f__Other',
            'p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae',
            'p__Proteobacteria;c__Other;o__Other;f__Other',
            'p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae',
            'p__Verrucomicrobiota;c__Other;o__Other;f__Other',
            'p__Other;c__Other;o__Other;f__Other',
            ],
        'large': 1,
        'large_groups': ['p__Firmicutes', 
                         'p__Bacteroidota', 
                         'p__Actinobacteriota',
                         'p__Proteobacteria', 
                         'p__Verrucomicrobiota'],
        'palette': {
            'p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae': '#7D3560',
            'p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae': '#A1527F',
            'p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Family_XI':  "#CC79A7",
            'p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae': "#E794C1",
            'p__Firmicutes;c__Other;o__Other;f__Other': "#EFB6D6",
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae': '#098BD9',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae': '#56B4E9',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae': '#7DCCFF',
            'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae': '#BCE1FF',
            'p__Bacteroidota;c__Other;o__Other;f__Other': '#E7F4FF',
            'p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae': '#6D9F06',
            'p__Actinobacteriota;c__Other;o__Other;f__Other': '#BDEC6F',
            'p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae': '#C17754',
            'p__Proteobacteria;c__Other;o__Other;f__Other': '#FCB076',
            'p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae': '#148F77',
            'p__Verrucomicrobiota;c__Other;o__Other;f__Other': '#48C9B0',
            'p__Other;c__Other;o__Other;f__Other': '#B7B7B7'
            },
    }
}

oral = {
    2: {'small'}
}


# ### FECAL BAR PLOT ###
fecal_tax_order_hard = [
    'p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae',
    'p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae',
    'p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Family_XI',
    'p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae',
    'p__Firmicutes;c__Other;o__Other;f__Other',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae',
    'p__Bacteroidota;c__Other;o__Other;f__Other',
    'p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae',
    'p__Actinobacteriota;c__Other;o__Other;f__Other',
    'p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae',
    'p__Proteobacteria;c__Other;o__Other;f__Other',
    'p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae',
    'p__Verrucomicrobiota;c__Other;o__Other;f__Other',
    'p__Other;c__Other;o__Other;f__Other'
]

feces_keep_phyla = ['p__Firmicutes', 'p__Bacteroidota', 'p__Actinobacteriota',
                    'p__Proteobacteria', 'p__Verrucomicrobiota']

feces_palette_hard = {
    'p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae': '#7D3560',
    'p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae': '#A1527F',
    'p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Family_XI':  "#CC79A7",
    'p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae': "#E794C1",
    'p__Firmicutes;c__Other;o__Other;f__Other': "#EFB6D6",
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae': '#098BD9',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae': '#56B4E9',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae': '#7DCCFF',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae': '#BCE1FF',
    'p__Bacteroidota;c__Other;o__Other;f__Other': '#E7F4FF',
    'p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae': '#6D9F06',
    'p__Actinobacteriota;c__Other;o__Other;f__Other': '#BDEC6F',
    'p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae': '#C17754',
    'p__Proteobacteria;c__Other;o__Other;f__Other': '#FCB076',
    'p__Verrucomicrobiota;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae': '#148F77',
    'p__Verrucomicrobiota;c__Other;o__Other;f__Other': '#48C9B0',
    'p__Other;c__Other;o__Other;f__Other': '#B7B7B7'
}


### Vaginal barplot ###
vaginal_genera = [
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus',
    'd__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella',
    'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Gardnerella',
    'd__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae;g__Megasphaera',
    'd__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Shuttleworthia',
    'd__Bacteria;p__Actinobacteriota;c__Coriobacteriia;o__Coriobacteriales;f__Atopobiaceae;g__Atopobium',
    'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium',
    'd__Other;p__Other;c_Other;o__Other;f__Other;g__Other',
]
vaginal_cmap = {
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus': '#66c2a5',
    'd__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella': '#fc8d62',
    'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Gardnerella': '#8da0cb',
    'd__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae;g__Megasphaera': '#e78ac3',
    'd__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Shuttleworthia': '#a6d854',
    'd__Bacteria;p__Actinobacteriota;c__Coriobacteriia;o__Coriobacteriales;f__Atopobiaceae;g__Atopobium': '#ffd92f',
    'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium': '#e5c494',
    'd__Other;p__Other;c_Other;o__Other;f__Other;g__Other': '#b3b3b3',
}

#### Nasal Barplots ###
nasal_colors = {
    'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Corynebacteriales;f__Corynebacteriaceae;g__Corynebacterium': '#6baed6',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus': '#a50f15',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Dolosigranulum': '#de2d26',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus': '#fb6a4a',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus': '#fcae91',
    'd__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae;g__Veillonella': '#fd8d3c',
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Moraxella': '#54278f',
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Pasteurellaceae;g__Haemophilus': '#756bb1',
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas': '#9e9ac8',
    'd__Other;p__Other;c_Other;o__Other;f__Other;g__Other': '#b3b3b3',
}

nasal_keep_class = ['c__Gammaproteobacteria', 'c__Bacilli', 'c__Negativicutes', 'c__Actinobacteria']

nasal_order = [
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Moraxella',
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Pasteurellaceae;g__Haemophilus',
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Dolosigranulum',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus',
    'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus',
    'd__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae;g__Veillonella',
    'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Corynebacteriales;f__Corynebacteriaceae;g__Corynebacterium',
    'd__Other;p__Other;c_Other;o__Other;f__Other;g__Other',
]

oral_order = [
    'p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae',
    'p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Gemellaceae',
    'p__Firmicutes;c__Bacilli;o__Lactobacillales;f__P5D1-392',
    'p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae',
    'p__Firmicutes;c__Bacilli;o__Other;f__Other',
    'p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Neisseriaceae',
    'p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Pasteurellaceae',
    'p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae',
    'p__Proteobacteria;c__Gammaproteobacteria;o__Other;f__Other',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae',
    'p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae',
    'p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Weeksellaceae',
    'p__Bacteroidota;c__Bacteroidia;o__Other;f__Other',
    'p__Firmicutes;c__Negativicutes;o__Veillonellales-Selenomonadales;f__Veillonellaceae',
    'p__Actinobacteriota;c__Actinobacteria;o__Micrococcales;f__Micrococcaceae',
    'p__Actinobacteriota;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae',
    'p__Actinobacteriota;c__Actinobacteria;o__Other;f__Other',
    'p__Other;c__Other;o__Other;f__Other',
    ]

oral_classes =['c__Bacilli', 
               'c__Gammaproteobacteria', 
               'c__Bacteroidia',
               'c__Negativicutes', 
               'c__Actinobacteria']