# QARM_RISK

Repository for the Quantitative asset and risk management - Risk Project

![HEC Lausanne Logo](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a3/HEC_Lausanne_logo.svg/293px-HEC_Lausanne_logo.svg.png)

This project is based on :


------------------------------------------------------------------------------------------------------------------------------
**DRAWDOWN: FROM PRACTICE TO THEORY AND BACK AGAIN**

Lisa R. Goldberg and Ola Mahmoud

 You can find this paper [here](https://arxiv.org/pdf/1404.7493.pdf).
 
------------------------------------------------------------------------------------------------------------------------------

This paper introduce the **Conditional Expected Drawdown** which is the tail mean of the distribution of Drawdown. This risk measure is coherent and therefore allows to construct portfolio based on it and to compute the marginal contribution of an asset to the overall risk. 

<img src="https://github.com/blacksouane/QARM_RISK/blob/master/Plots/MDD_Distribution.png" height="400" width="600"> 

Hence, the purpose of this project is to expose the properties of maximum Drawdown and CED. Typically, we will show the relation between CED and the frequency of observation, the auto-correlation of the underlying process and the length of the path.

We will then look at the relation between CED, the Drawdown, the speed of the drawdown and the recovery period.

Finally, we use different asset allocation methods and compute the contribution to CED of each asset. The allocation are : 

1. Equally Weighted
2. Risk Parity

The last part of the paper takes intraday bitcoin prices to look at the properties of CED for high frequency trading and look at the properties of a CED parity allocation. 

## Data

We use two dataset of the Fama/French data:

- Fama/French [3 factors](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/f-f_factors.html) daily (to extract the rf and rm values)
- Fama/French [10 industry Portfolio](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html) daily

We also used intraday bictoin prices. The data are to heavy to be on this repository but can find them here : [Bitcoin Data on Kaggle](https://www.kaggle.com/mczielinski/bitcoin-historical-data). The code should work even if you time span is bigger than ours.

## Implementation

The code is relatively computotionnaly intensive and therefore takes a lot of time to run (around 10 hours on my computer). Thus, we have provided a instance of our Workspace with all the data and computations already done. This can allow you to just do or change to computations that you want and not everything.

### Authors

* **Antoine-Michel Alexeev** 
* **Julien Bisch**
* **Benjamin Souane** 
* **Ludovic Suchet**

### Acknowledgement

* **Prof. Eric Jondeau and Alexandre Pauli** for the help provided.
* Thanks to anyone whose code was used.
