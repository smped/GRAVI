# To-Do List

## Input Validity Functions

## Motifs

- Include motif enrichment as part of the standard workflow
- Include pathway enrichment

The problematic range is `chr11:60659881-60660540`.
This appears mapped to different features in both AR & ER:

```
$AR
GRanges object with 1 range and 17 metadata columns:
      seqnames            ranges strand | n_windows      n_up    n_down            keyval_range   AveExpr     logFC     P.Value
         <Rle>         <IRanges>  <Rle> | <integer> <integer> <integer>               <GRanges> <numeric> <numeric>   <numeric>
  [1]    chr11 60660001-60660450      * |         7         2         0 chr11:60660101-60660250 -0.733593   1.54097 8.69699e-06
              fdr overlaps_ref direction                                             gene_id             gene_name regulatory_feature
        <numeric>    <logical>  <factor>                                     <CharacterList>       <CharacterList>           <factor>
  [1] 2.23525e-05        FALSE        Up ENSG00000172689,ENSG00000110104,ENSG00000149506,... MS4A10,CCDC86,ZP1,...     Super Enhancer
      gene_region ihw_covariate     fdr_ihw   status
         <factor>      <factor>   <numeric> <factor>
  [1]   Gene Body  ER+ H3K27ac+ 2.23766e-05       Up
  -------
  seqinfo: 24 sequences from GRCh37 genome

$ER
GRanges object with 1 range and 17 metadata columns:
      seqnames            ranges strand | n_windows      n_up    n_down            keyval_range   AveExpr     logFC   P.Value       fdr
         <Rle>         <IRanges>  <Rle> | <integer> <integer> <integer>               <GRanges> <numeric> <numeric> <numeric> <numeric>
  [1]    chr11 60659881-60660540      * |         9         1         0 chr11:60660061-60660240 0.0403986   0.65486 0.0958111  0.135919
      overlaps_ref direction                                             gene_id                 gene_name regulatory_feature
         <logical>  <factor>                                     <CharacterList>           <CharacterList>           <factor>
  [1]         TRUE        Up ENSG00000110104,ENSG00000149506,ENSG00000257052,... CCDC86,ZP1,AP003721.4,...           Enhancer
      gene_region ihw_covariate   fdr_ihw    status
         <factor>      <factor> <numeric>  <factor>
  [1]   Gene Body      H3K27ac+  0.127014 Unchanged
  -------
  seqinfo: 24 sequences from GRCh37 genome
```

The range ER maps to is immediately upstream and doesn't overlap!

```
> feat_gr %>% subset(vapply(gene_name, function(x){"CCDC86" %in% x}, logical(1)))
GRanges object with 15 ranges and 3 metadata columns:
       seqnames            ranges strand |     feature                                             gene_id                  gene_name
          <Rle>         <IRanges>  <Rle> | <character>                                     <CharacterList>            <CharacterList>
   [1]    chr11 60507308-60513467      * |    Enhancer ENSG00000181995,ENSG00000166959,ENSG00000214782,... LINC00301,MS4A8,MS4A18,...
   [2]    chr11 60518561-60519157      * |    Enhancer ENSG00000181995,ENSG00000166959,ENSG00000214782,... LINC00301,MS4A8,MS4A18,...
   [3]    chr11 60588922-60590194      * |    Enhancer ENSG00000214782,ENSG00000166961,ENSG00000172689,...   MS4A18,MS4A15,MS4A10,...
   [4]    chr11 60593031-60593270      * |    Enhancer ENSG00000214782,ENSG00000166961,ENSG00000172689,...   MS4A18,MS4A15,MS4A10,...
   [5]    chr11 60593620-60594148      * |    Enhancer ENSG00000214782,ENSG00000166961,ENSG00000172689,...   MS4A18,MS4A15,MS4A10,...
   ...      ...               ...    ... .         ...                                                 ...                        ...
  [11]    chr11 60639083-60639504      * |    Enhancer ENSG00000172689,ENSG00000110104,ENSG00000149506,...      MS4A10,CCDC86,ZP1,...
  [12]    chr11 60641711-60642400      * |    Enhancer ENSG00000172689,ENSG00000110104,ENSG00000149506,...      MS4A10,CCDC86,ZP1,...
  [13]    chr11 60655203-60662802      * |          SE ENSG00000172689,ENSG00000110104,ENSG00000149506,...      MS4A10,CCDC86,ZP1,...
  [14]    chr11 60684222-60684828      * |    Enhancer ENSG00000110104,ENSG00000149506,ENSG00000257052,...  CCDC86,ZP1,AP003721.4,...
  [15]    chr11 60707129-60707728      * |    Enhancer ENSG00000110104,ENSG00000149506,ENSG00000257052,...  CCDC86,ZP1,AP003721.4,...
  -------
  seqinfo: 24 sequences from GRCh37 genome
```
