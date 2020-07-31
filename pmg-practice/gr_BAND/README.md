# Example of band structure pymatgen parsing on graphene

* Run VASP self-consistently with relaxation (IBRION = 2, NSW = N>1). I did NSW = 80.
* Run VASP non self-consistently without relaxation and with fixed charge density (ICHARG = 11, IBRION = -1, NSW = 1).
* Check analysis.ipynb for parsing and analysis.