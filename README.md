
A Python script for 2 dimensional grain packing and permeability calculations is provided as a reference. For a detailed explanation of the algorithm, refer to the study "Relationship between Grain Packing and Porosity," which is available at https://arxiv.org/pdf/ and was also collected in the Catalogue MSC AT PORELAB 2019 at https://porelab.no/wp-content/uploads/2020/04/Catalogue_MSc_final_oppslag.pdf.

<img src="CatalogueMasterStudentsPorelab2019Opportunities2020_QiangqiangMao.jpg" alt="BriefIntroduction" style="width: 75%;" align="center">

This code was made for the course TPG4560 - Petroleum Engineering, Specialization Project in Norwegian University of Science and Technology(NTNU) as an attenmpt to understand the relationship between grain packing and single phase transport properties. 

In this repository, the codes for generating 2 dimensional grain packing, area weighted normal distribution and permeability calculation based on Poisson equation are presented on the python platform. Specifically, the awcdf.py is the code for area weighted normal distribution to determine grain size distribution. the GrainPacking_2DMain.py is the main code for generating 2D grain packing. the GrainPacking2D_Function.py is the code for used function that is assistance of the GrainPacking_2DMain.py. 

for 2D grian packing permeability, we mainly calculate it by solving Poisson equation using Finite difference Method. The code is given in permeabilityFunction.py. And the mesh.npy is one example of the packing data which will be used to calculate permeability. 

The code will run fastly. The grain packing process will take around half minite and permeability calculation will take about 5 seconds.

Bear in mind that the code is far from perfect. Some code optimization here and there can be done. Especially, there are probably some bugs I haven't found. If someone find bugs, do not hesitate to give comment and let me know. It will be appreciated that you are willing to make some improvements. Any suggestion and modification are warmly welcome. you can send it to my email, upcmaoqiangqiang@outlook.com or maoqiangqiangupc@gmail.com. Thank you.

Sincerely thanks are given to my supervisor Prof. Carl Fredirk Berg in the department of geoscience and petroleum in NTNU. He gave me many guidance on my graduation specialized project and my bachelor thesis. Thanks my supervisor.

by Mao Qiangqiang 4th, July 2019, exchanged student in Norway University of Science and Technology and bachelor in China University of Petroleum (East China).
