# QARM_RISK

Repository for the Quantitative asset and risk management - Risk Project

![HEC Lausanne Logo](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a3/HEC_Lausanne_logo.svg/293px-HEC_Lausanne_logo.svg.png)

This project is based on :


* **DRAWDOWN: FROM PRACTICE TO THEORY AND BACK AGAIN**
Lisa R. Goldberg and Ola Mahmoud
You can find this paper [here](https://arxiv.org/pdf/1404.7493.pdf).

This paper introduce the **Conditional Expected Drawdown** which is the tail mean of the distribution of Drawdown. This risk measure is coherent and therefore allows to construct portfolio based on it and to compute the marginal contribution of an asset to the overall risk. 

![Schema_CED](https://github.com/blacksouane/QARM_RISK/blob/master/Plots/MDD_Distribution.png)

Hence, the purpose of this project is to expose the properties of maximum Drawdown and CED. Typically, we will show the relation between CED and the frequency of observation, the auto-correlation of the underlying process and the length of the path.

## Implementation

The code is relatively computotionnaly intensive and therefore takes a lot of time to run (around 10 hours on my computer). Thus, we have provided a instance of our Workspace with all the data and computations already done. This can allow you to just do or change to computations that you want and not everything.

### Authors

* **Alessia Di Pietro** 
* **Antoine-Michel Alexeev** 
* **Benjamin Souane** 

### Acknowledgement

* **Prof. Eric Jondeau and Alexandre Pauli** for the help provided.
* Thanks to anyone whose code was used.
