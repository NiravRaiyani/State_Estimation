# State_Estimation
Implementation of Kalman Filter, Extended Kalman Filter and Moving Horizon Estimation to the stirred tank mixing process. This repository uses the same system as the one used in [Implementation and comparison of Advanced process control to stirred tank mixing process](https://github.com/NiravRaiyani/State_Estimation).


## Stirred Tank Mixing Process
![](assets/mixing.png)

The Continuously stirred mixing tank in the figure has two time-varying inlets F1(t) and F2(t) with
different density. The density of both the inlets is constant and is given by rho1 and rho2. It is assumed
that the tank is well mixed so that outlet F(t) has the same density as the density in the tank i.e.
rho(t). The volume of the tank occupied by the liquid is V(t) and the corresponding height is h(t)
with the surface area S. 
## Folder Contents 
* KF

```KF_3.m```Main File. 

* EKF:

```EKF1.m``` Main file.

```lin1.m``` To compute Jacobian at different states.

* MHE

```MHE_3.m```Main File.

```MHE_cost```To calculate the total cost at every time step.

```lin1.m```To compute Jacobian at different states.

```mhe.m```To perform constrained non-linear optimization for state estimation.

```sys.m```To calculate one step ahed prediction.

Code flow for MHE:
![](assets/pdf1.png)
