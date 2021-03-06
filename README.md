# Relative Sink-Strenght

## Goal:
- Compute partitioning factors with the relative sink-strength principle.

## Scientific Bases:

- https://academic.oup.com/jxb/article/61/8/2203/488729
- https://academic.oup.com/aob/article/91/3/361/231163

## Example:
- At this repo you will find an example ([relative_ss_example.R](relative_ss_example.R)) with an arbitrary crop composed by 3 organs: 
  - 1 root
  - 1 leaf
  - 1 stem
- We also hypothesize a scenario where water stress can alter the crop's sink strength where each organs has it's own sensitivity. If we consider that the roots' SS is more resilient to water stress then it'd mimic a higher partitioning factor to roots when the crop is subjected to drought.

## Partitioning Factors under optimum water
![alt text](https://github.com/Murilodsv/relative_ss/blob/master/pfac_optimum_h2o.png)

## Partitioning Factors under water Stress
![alt text](https://github.com/Murilodsv/relative_ss/blob/master/pfac_h2o_stress.png)
