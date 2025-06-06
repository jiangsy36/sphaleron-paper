# sphaleron-paper

In this work, I reproduce the main results of the paper arXiv:2505.05608 and arXiv:2307.04713.

In the folder "FDM", I solve the time evolution differential equations of sphaleron by using the fourth order Runge-Kutta method. The results are consistent with the Fig.3 of 2505.05608.
![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/FDM/fA.png)
![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/FDM/fB.png)
![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/FDM/fC.png)
![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/FDM/K.png)
![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/FDM/H.png)
In the folder spectral, we solve the sphaleron by using the spectral method where the solutions are expanded by Chebyshev Polynomials. The caculation details are given in the appendix of 2307.04713. The sphaleron differential equations are given by "A.19-A.21" (after axis transformation):


![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/spectral/sphaleron.png)
The solution reads,
![image](https://github.com/jiangsy36/sphaleron-paper/blob/main/spectral/hdelta.png)
