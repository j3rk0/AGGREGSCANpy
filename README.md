# AGGREGSCANpy
AGGREGSCAN is a script to compute protein aggregation profile features.

this repository is a python implementation of:
"Conchillo-Solé, O., de Groot, N.S., Avilés, F.X. et al. AGGRESCAN: a server for the prediction and evaluation of "hot spots" of aggregation in polypeptides. BMC Bioinformatics 8, 65 (2007). https://doi.org/10.1186/1471-2105-8-65" 

# Example usage

```python
from aggregscan import aggregation_profile

sequence = "MMCAATSPGKHQAATSPGKHQNRDEIVFGKHQNRDE"  
result = aggregation_profile(sequence)
print(result)
```

the output of aggregation profile function is a dictionary. for the previous example the output should be:

```python
{'HST': -0.02, 
'a3vSA': -0.46927777777777785, 
'nHS': 1, 
'NnHS': 2.7777777777777777, 
'AAT': 3.0813, 
'THSA': 0.8876000000000002, 
'TA': -16.548399999999997, 
'AAATr': 0.08559166666666668, 
'THSAr': 0.02465555555555556, 
'Na4vSS': -49.61, 
'Hot Spots': [(0, 4)], 
'HSAs': [0.8876000000000002], 
'A4VHS': [0.25805000000000006], 
'a3v': [ 0.91 ,  0.91 ,  0.604, -0.036, -0.036, -0.159, -0.294, -0.334,
       -0.535, -0.931, -1.033, -1.231, -0.036, -0.036, -0.159, -0.294,
       -0.334, -0.535, -0.931, -1.033, -1.231, -1.302, -1.24 , -1.836,
       -1.412,  1.822,  1.594,  1.754, -0.535, -0.931, -1.033, -1.231,
       -1.302, -1.24 , -1.836, -1.412], 
'a4v': [ 0.1526,  0.1526,  0.4704,  0.2566,  0.0158, -0.1718, -0.2716,
       -0.4506, -0.6254, -0.8128, -0.7532, -0.6534, -0.499 , -0.3512,
       -0.1718, -0.2716, -0.4506, -0.6254, -0.8128, -1.0064, -1.1474,
       -1.3284, -1.4042, -0.7936, -0.2144,  0.3844,  0.6446,  0.7408,
        0.1698, -0.3952, -1.0064, -1.1474, -1.3284, -1.4042, -1.375 ,
       -1.375 ]}
```

feel free to check output consistency with http://bioinf.uab.es/aggrescan/ .
For any information about the method and the features interpretation refer to the original paper.
