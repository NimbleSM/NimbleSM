# J2 plasticity model in NimbleSM

## Brief description

The model is an extension of the classical J2 plasticity for small strains to the finite deformation regime. It does this by means of the logarithmic and exponential maps, and thus its elasticity is of the Henky type. It includes power law hardening and power law rate dependence (viscoplasticity).

## Material parameters

A typical line for J2 plasticity with its input parameters is:
```
material parameters:              material_1 j2-plasticity elastic_modulus 130.6e+09 poissons_ratio 0.343 density 8960 yield_stress 89.7e+06 hardening_exponent 0 reference_plastic_strain 0.02 reference_viscoplastic_stress 89.7e+06 rate_dependence_exponent 0 reference_plastic_strain_rate 1e+04
```
where
* `elastic_modulus`
  * Elastic or Young's modulus
* `poissons_ratio`
  * Poisson's ratio
* `density`
  * Reference mass density
* `yield_stress`
  * Reference yield of flow stress
* `hardening_exponent`
  * Exponent for power law hardening. If zero there is no hardening.
* `reference_plastic_strain`
  * Reference plastic strain for power law hardening
* `reference_viscoplastic_stress`
  * Reference viscoplastic stress for power law rate dependence. Typical value is the same as `yield_stress`. If zero there is no rate dependence.
* `rate_dependence_exponent`
  * Exponent for power law rate dependence. If zero there is no rate dependence.
* `reference_plastic_strain_rate`
  * Reference plastic strain rate for power law rate dependence.