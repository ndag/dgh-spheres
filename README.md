# dgh-spheres


This GitHub repository contains code complementing the paper 

``The Gromov-Hausdorff distance between spheres", by Sunhyuk Lim, Facundo Mémoli and Zane Smith. (to appear, Geometry and Topology. The [arxiv version is here](https://arxiv.org/abs/2105.00611).)


<img src="/dm-spheres-mat.png" width="500" height="auto">

## The correspondence from Appendix D of the [arxiv version](https://arxiv.org/abs/2105.00611))
The matlab script *corresp_s1_s2_app_D.m* provides code for constructing the correspondence between S1 and S2 which is describdin the Appendix D of the ArXiv version of the paper. The code also computes its distortion.



# Embedding Projection Corresponds for estitating the Gromov-Hausdorff distance between spheres

This code implements certain surjective maps between spheres which have 'low' distortion. As such these can be used to upper bound the **Gromov-Hausdorff distance** between spheres (when endowed with their geodesic distances).

These ideas were first proposed by **Zane Smith**. In particular, circa 2015/2016 Zane proposed the surjective maps $S^1 \to S^2$, $S^1 \to S^3$, $S^1 \to S^4$ and $S^1 \to S^5$ and extensively tested their distortion experimentally. 

See [3] for a historical account.

The critical idea that Zane had was to use the closest point projecton of a given sphere onto the image of an embedding of the lower dimensional sphere into the larger dimensional sphere. 

For example, in the case of $S^1$ versus $S^2$, Zane's idea was to first consider the curve $\tilde{\gamma}_2(t) = (\cos(t),\sin(t),\cos(3t))$, for $t\in[0,2\pi]\cong S^1$ and then define $\gamma(t) = \frac{\tilde{\gamma}_2(t)}{\|\\gamma_0(t)\|}$ as curve lying on the surface of $S^2$. Note that $\gamma:S^1 \to S^2$ is a diffeo onto its image.

Then, if $P:S^2\rightarrow \mathrm{trace}(\gamma)$ is the closest point projection map, Zane deefined the surjection
$\Psi_{2,1}:S^2\rightarrow S^1, \mbox{s.t. $x \mapsto \gamma^{-1}(P(x))$}.$

Zane experimentally verified that the distortion of $\Psi_{2,1}$ is $\frac{2\pi}{3}$.

In a similar way, by instead considering a curve the curve $\gamma_3(t) = (\cos(t),\sin(t),\cos(3t),\sin(3t))$, Zane constructed a surjection $\Psi_{3,1}:S^3 \to S^1$ which he studied experimentally to conclude that its distortion is $\frac{2\pi}{3}$ as well.


A generalization to the case of $S^1\to S^n$ is as follows:

when $n=2k+1$, $$\tilde{\gamma}_n(t) = (\cos(t),\sin(t),\cos(3t),\sin(3t),\ldots,\cos((2k+1)t),\sin((2k+1)t));$$

when $n=2k$,   $$\tilde{\gamma}_n(t) = (\cos(t),\sin(t),\cos(3t),\sin(3t),\ldots,\cos((2k+1)t)).$$



## For background, see:

[1] The Gromov-Hausdorff distance between spheres. Sunhyuk Lim, Facundo Memoli and Zane Smith. Geometry and Topology. To appear. https://arxiv.org/abs/2105.00611

[2] Gromov-Hausdorff distances, Borsuk-Ulam theorems, and Vietoris-Rips complexes. Henry Adams, Johnathan Bush, Nate Clause, Florian Frick, Mario Gómez, Michael Harrison, R. Amzi Jeffs, Evgeniya Lagoda, Sunhyuk Lim, Facundo Mémoli, Michael Moy, Nikola Sadovek, Matt Superdock, Daniel Vargas, Qingsong Wang, Ling Zhou. https://arxiv.org/abs/2301.00246

[3] Embedding-Projection Correspondences for the estimation of the Gromov-Hausdorff distance.
F. M\'emoli and Z. Smith. 
