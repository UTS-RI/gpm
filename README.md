Author: le.gentil.cedric@gmail.com

## Gaussian Preintegrated Measurements (GPMs)

This repository contains a Matlab implementation of our paper __Gaussian Process Preintegration for Inertial-Aided State Estimation__ \[[Video](https://www.youtube.com/watch?v=T5rXOuAszmc)\].
It is a novel continuous and probabilistic preintegration method that leverages Gaussian Process regression and the application of linear operators to the covariance functions.


__ATTENTION__: Dear Mac users, it seems that there is an "issue" with Matlab and Catalina (and maybe other versions) as the hyper-parameter training from the GPML Toolbox can diverge (therefore leading to erroneous GPMs). In such case, please change `ncg =50` on line 22 of `Gpm.m` to `1` or `2` (that limits the number of iterations of the GPML hyper-parameter training), or simply use the hyper-parameter initial guess by setting `max_nb_loops` to `0` on  line 358 of `Gpm.m`.
The results would be suboptimal, but the overall method should work. Do not hesitate to contact me if the error persists.



#### Disclaimer

This code has been developed in a research context. It is not optimised in any way, and bugs might occur. The naming convention does not match the standard convention used for Matlab scripting. The current documentation is fairly light regarding the functions input/output specifications (the variables' name should be enough in most cases).
Nonetheless, I hope that this code, released under the [MIT License](LICENSE), can be useful to someone.
I will welcome any feedback and consider modifications/improvements as much as my schedule let me to.


#### Usage

This code relies on the [GPML Toolbox](http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html) developed by C.E. Rasmussen and H. Nickisch.
The present repository includes this toolbox.
Additionally, this code requires the Robotics System Toolbox (for rotation conversion functions).
Let me know if other dependencies should be mentioned here.

The function that actually computes the __GPMs__ is present in `Gpm.m`.
A basic demonstration script is contained in `gpm_demo.m`.
Here is an example of output:
```
PM error: 0.63825 m    1.3545 deg
GPM error: 0.050035 m    0.11747 deg
```

The files `gpm_jacobian_test.m` and `gpm_variance.m` provide tools to compare both the GPMs' postintegration correction Jacobians, and the GPMs' variances, with the expected value computed numerically.

The file `dev.m` shows the code used to derive the various kernel integrals and differentiation.

The other files contain auxiliary code needed in the scripts/function mentioned above.

The type of motion simulated (the aggressiveness _'slow'_ / _'fast'_, and rotations around 1-axis( _true_ )/multi-axis( _false_ )) can be changed by modifying the `traj_profile` and `one_axis` attributes of the simulation option argument.

#### Implementation details and known issues

###### Kernel choice
In this implementation, I used a _square exponential_ kernel function as it is simple to derive, and its hyper-parameters have intuitive meanings. Other kernels can easily be used instead.

###### 1-axis rotation case
Both the 1-axis and the multi-axis rotations cases have been implemented. The 1-axis formulation in our paper allows the analytical integration of the angular velocities when the system's orientation is confined to one fix rotation axis. In real scenarios, the IMU is subject to additive biases. With a 1-DoF gyroscope, the bias does not perturb the integration. However, with a 3-DoF gyroscope, the presence of unknown biases can break the assumption of rotation commutativity needed to analytically integrate the 3-axis angular velocities. Consequently, when dealing with 3-axis gyroscopes (in the case of 1-axis rotations), one needs to know the axis of rotation. Here (in the case of 1-axis rotations) we consider the rotations are around the z-axis.

###### Rotation covariance in multi-axis rotation case
In the multi-axis rotation case, the iterative noise propagation from [C. Forster et al.](http://rpg.ifi.uzh.ch/docs/TRO16_forster.pdf) is used to compute the covariance of the rotational part of the GPM.
I am not entirely sure why (implementation bug not excluded), but the variances estimated this way are one or two orders of magnitude too ...certain.... --> _To be investigated._


##### TODO

- Document the functions input/output
- Investigate the rotational variance in the multi-axis rotation scenario



#### Citing
This work has been published in the IEEE Robotics and Automation Letters and presented at ICRA 2020 Paris.
Please cite the following paper if you are using GPMs in your work:

> C. Le Gentil, T. Vidal-Calleja, and S. Huang, “__Gaussian Process Preintegration for Inertial-Aided State Estimation__,” IEEE Robotics and Automation Letters, vol. 5, no. 2, pp. 2108–2114, 2020.

```bibtex
@article{LeGentil2020,
	title = {{Gaussian Process Preintegration for Inertial-Aided State Estimation}},
	author = {{Le Gentil}, Cedric and Vidal-calleja, Teresa and Huang, Shoudong},
	journal = {IEEE Robotics and Automation Letters},
	number = {2},
	pages = {2108--2114},
	volume = {5},
	year = {2020} 
}


```

You will find the paper in the [IEEE Xplore](https://ieeexplore.ieee.org/document/8979155) and a video on [Youtube](https://www.youtube.com/watch?v=T5rXOuAszmc).
