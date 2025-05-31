# Generalized nonlinear Schrödinger equation

## 📜 **Generalized NLSE Used in the Code**

The evolution of the complex envelope $A(x, y, t; z)$ along the propagation axis $z$ is governed by:

$$
\boxed{
\begin{aligned}
\frac{\partial A}{\partial z} &= \underbrace{i \frac{1}{2k_0} \left( \frac{\partial^2 A}{\partial x^2} + \frac{\partial^2 A}{\partial y^2} \right)}_{\textbf{(1) Diffraction}} 
\underbrace{-i \frac{\beta_2}{2} \frac{\partial^2 A}{\partial t^2} + \frac{\beta_3}{6} \frac{\partial^3 A}{\partial t^3}}_{\textbf{(2) Dispersion (GDD + TOD)}}
 \underbrace{-i \gamma \left(1 + \frac{i}{\omega_0} \frac{\partial}{\partial t}\right)\left(|A|^2 A\right)}_{\text{(3) Kerr + self-steepening}} 
 \underbrace{-\frac{\alpha}{2} A}_{\text{(4) Linear absorbtion}} 
\end{aligned}
}
$$

where:

- $ A(x,y,t;z) $: complex field envelope
- $ k_0 = \omega_0 / c $: central wavenumber
- $ c $: speed of light in vacuum
- $ \beta_2, \beta_3 $: group delay dispersion (GDD) and third-order dispersion (TOD)
- $ \gamma $: Kerr nonlinear coefficient
- $ \omega_0 $: central angular frequency
- $ \alpha $: linear absorption coefficient

This makes our system capable of modeling ultrashort laser pulse propagation in **nonlinear, lossy media**. There are two ways to derive the GNLSE: from the Agrawal book and from the Ursula Keller book, because it depends on the Fourier transform convention. I choosed Ursula Keller's notation, because python's FFT and IFFT library incorporates the Ursula Keller's convention.

---