# PanIE
**Pan**genome **I**ndependent **E**dge identification

The repository contains python scripts to create a non-redundant set of graph edges to represent all variants in a pangenome graph. 

## Overview
In pangenome graphs, a single genetic variant is often represented by multiple edges. Treating all of these as independent signals inflates multiple testing and give redundant edges extra chances to ‘win’ an association, which bias downstream association analyses (e.g., GWAS/eQTL).

To remove the redundancy, our method identifies a single best representation edge for each variant instead of picking one randomly for each variant. In our study, we prioritize edges with higher read depths (although users can adjust the prioritization criteria as needed). This approach accounts for the fact that when multiple edges represent the same variant, they may not represent the variant equally well due to differences in sequence complexity.

Our edge filtering method effectively reduces the multiple testing burden, and at the same time, retains high quality edges that are most informative for downstream eQTL/GWAS analysis, which is an important step to ensure pangenome specific eQTL discovery.
<p align="center">
<img width="492" height="148" alt="Screenshot 2025-11-21 at 1 13 48 AM" src="https://github.com/user-attachments/assets/3321d132-94b7-4c96-8b5a-24939c940080" />
</p>

## Approach
We decompose the pangnoeme graph into biconnected components to parallel the task and identify independent edges within each subgraph. Then, we iteratively select the best representative edge, prioritizing by higher read depth, and designate it as the representative edge for that variant. After selection, we remove its dependent edges from the graph. Dependent edges are defined as bridges in the graph, whose removal increases the number of connected components.
<p align="center">
<img width="503" height="173" alt="Screenshot 2025-11-21 at 1 08 51 AM" src="https://github.com/user-attachments/assets/8ffaf117-6e29-405d-80f6-5ec6bfe4c32f" />
</p>


## Quick Start
```bash
# 1. step1 - decompose pangenome graph into biconnected components
for chr in chr{1..22} chrX; do
    python filter_redundancy_bubble_puncture.fix2.step1_subgraph.py $chr 
done 
# 2. step2 - identify independent edge in each biconnected subgraph
python filter_redundancy_bubble_puncture.fix2.step2_filter.py $chr <subgraph.gml>

# 3. step3 (optional) - speed up independent edge identification by only check bridge edge within 50 steps 
# a patch/compromise in case that a subgraph is too complex and it takes too long to run step2
python filter_redundancy_bubble_puncture.fix2.step2_filter.dp50.py $chr <subgraph.gml>

```
### input 

### output

