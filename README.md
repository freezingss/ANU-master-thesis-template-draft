- Finished the full Laplace-EM derivation for FA-DMR, including:
  the design decision to run the E-step directly on the Q-dimensional group random effect rather than the K-dimensional factor space,
  the latter leads to a mathematically incorrect group-collapsing approximation and a systematic inflation bias in B via Laplace shrinkage of the posterior second moment.

- Main technical contribution: a Woodbury acceleration built on the Poisson-Multinomial equivalence. Introducing the auxiliary offset ξ decouples the categories in the E-step Hessian, bringing the per-Newton-step cost down from O(Q³) to O(QK²+K³),


- Proved that using the Poisson surrogate for the Newton direction, while checking the true multinomial log-posterior in the backtracking line search, still yields exact linear convergence to the true posterior mode.

- Separately characterized the systematic gap between the surrogate and true Hessians at that mode:
  the surrogate understates posterior uncertainty exactly in the directions controlled by σ²,
  motivated a hybrid scheme using the cheap surrogate for the search, then recompute the true curvature once at the mode.
