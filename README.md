# sphaleron-paper

In this work, I reproduce the main results of the paper arXiv:2505.05608 and arXiv:2307.04713.

In the folder "FDM", I solve the time evolution differential equations of sphaleron by using the fourth order Runge-Kutta method. The results are consistent with the Fig.3 of 2505.05608.

In the folder spectral, we solve the sphaleron by using the spectral method where the solutions are expanded by Chebyshev Polynomials. The caculation details are given in the appendix of 2307.04713. The sphaleron differential equations are given by "A.19-A.21":

$$
\begin{align}
	2\left(1+x\right)_{}^2 \bar{f}'' &= \left(2\bar{f}+1+x\right)\left(2\bar{f}-1+x\right)\left(2\bar{f}+x\right)\nonumber\\
	&+\frac{a_{}^2}{16\beta}\left(1+x\right)_{}^2\left(2\bar{f}-1+x\right)\left(2\bar{h}+1+x\right)_{}^2\nonumber\\
	&+\frac{a_{}^2 \varrho_{3}}{6\beta}\left(1+x\right)_{}^2\left(2\bar{f}-1+x\right)\left(2\bar{h}^{}_\Delta+1+x\right)_{}^2\;,\\
	\label{eq:hbartypeII}
	\left(1+x\right)_{}^2 \bar{h}'' + \left(1+x\right)\left(2\bar{h}'+1\right)&= \frac{1}{4}\left(2\bar{f}-1+x\right)_{}^2\left(2\bar{h}+1+x\right)\nonumber\\
	&-\frac{a_{}^2}{8\beta}\left(1+x\right)_{}^2\left\{\left(\varrho^{}_{1}-\varrho^{}_{2}\right)\left(2\bar{h}+1+x\right) \left[4-\left(2\bar{h}+1+x\right)_{}^2\right]\right.\nonumber\\
	&\left.-\varrho^{}_{2}\left(2\bar{h}+1+x\right)\left[\left(2\bar{h}+1+x\right)_{}^2-2\left(2\bar{h}^{}_\Delta+1+x\right)\right]\right.\nonumber\\
	&\left.+4\left(\varrho^{}_{4}-\varrho^{}_{1}+\varrho^{}_{2}\right)\left(2\bar{h}+1+x\right)\right.\nonumber\\
	&\left.+2\left(\sqrt{2\varrho^{}_{2}\varrho^{}_{3}\varrho^{}_{5}}-\varrho^{}_{2}\right)\left(2\bar{h}+1+x\right)\left(2\bar{h}^{}_\Delta+1+x\right)\right.\nonumber\\
	&\left.+\left(\varrho^{}_{1}-\varrho^{}_{4}-\sqrt{2\varrho^{}_{2}\varrho^{}_{3}\varrho^{}_{5}}\right)\left(2\bar{h}+1+x\right)\left(2\bar{h}^{}_\Delta+1+x\right)_{}^2\right\}\;,\\
	\label{eq:hbarDeltatypeII}
	\varrho^{}_{3} \left(1+x\right)_{}^2 \bar{h}_\Delta'' + \varrho^{}_{3} \left(1+x\right)\left(2\bar{h}_\Delta'+1\right)	&= \frac{2\varrho^{}_{3} }{3}\left(2\bar{f}-1+x\right)_{}^2\left(2\bar{h}^{}_\Delta+1+x\right)\nonumber\\
	&-\frac{a_{}^2 \varrho^{}_{2}}{8\beta} \left(1+x\right)_{}^2\left[\left(2\bar{h}+1+x\right)_{}^2-2\left(2\bar{h}^{}_\Delta+1+x\right)\right]\nonumber\\
	&+\frac{a_{}^2}{8\beta}\left(1+x\right)_{}^2\left[2\left(2\varrho^{}_{3}\varrho^{}_{5}-\varrho^{}_{2}\right)\left(2\bar{h}^{}_\Delta+1+x\right)\right.\nonumber\\
	&\left.-\left(\sqrt{2\varrho^{}_{2}\varrho^{}_{3}\varrho^{}_{5}}-\varrho^{}_{2}\right)\left(2\bar{h}+1+x\right)^2\right.\nonumber\\
	&\left.-\left(\varrho^{}_{1}-\varrho^{}_{4}-\sqrt{2\varrho^{}_{2}\varrho^{}_{3}\varrho^{}_{5}}\right)\left(2\bar{h}+1+x\right)_{}^2
	\left(2\bar{h}^{}_\Delta+1+x\right)\right.\nonumber\\
	&\left.+\left(\varrho^{}_{1}-\varrho^{}_{4}-\varrho^{}_{3}\varrho^{}_{5}-\sqrt{\varrho^{}_{2}\varrho^{}_{3}\varrho^{}_{5}/2}\right)\left(2\bar{h}^{}_\Delta+1+x\right)_{}^3\right]\;.
\end{align}
$$
