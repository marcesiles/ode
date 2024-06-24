# Welcome to ODE

In this project, we delve into the exploration of three distinct numerical methods for solving differential equations. The methods we will examine are:

1. **Euler Method**: A straightforward and fundamental approach, providing a basic understanding of numerical integration techniques.
2. **Runge-Kutta Method of Order 2**: An intermediate method that offers improved accuracy over the Euler method by considering additional points within the interval of interest.
3. **Runge-Kutta Method of Order 4**: Highly accurate and widely used method, balancing computational efficiency and precision by averaging the slopes at multiple points of interest.

Through this exploration, we aim to undestand the principles, implementation and comparative effectiveness of these numerical methods when it comes to solve ordinary differential equations.

## Euler Method

The Euler method is a straightforward approach based on the Taylor expansion of the function $x(t)$. The expansion is expressed as:

$$
\text{Taylor Expansion} \Rightarrow x(t+h) = x(t) + h\frac{dx}{dt} + \overbrace{\frac{h^2}{2} \frac{d^2x}{dt^2}}^{\varepsilon} + \mathcal{O}(h^3).
$$

To advance in time by a step $h$, which we assume to be sufficiently small, we use the equation:

$$
\boxed{x(t + h) = x(t) + h\cdot f(x,t).}
$$

The error associated with this approximation is related to the number of time steps $N$. This associated error can be estimated as follows:

\begin{align}
\sum\varepsilon = \sum_{k=0}^{N-1}\frac{h^2}{2}\left. \frac{d^2x}{dt^2} \right|{x_{k}, t_{k}} = \frac{h}{2}\sum{k=0}^{N-1}h\left.\frac{df}{dt}\right|{x_{k}, t_{k}}
\\
\approx \frac{h}2\int_{a}^b\frac{df}{dt}d t = \frac{h}{2}\left[f_{b} - f_{a}\right].
\end{align}

In this equation, we assume taking $N = (b-a)/h$ time steps to reach the final point. Therefore, the total approximation error depends linearly on $h$ multiplied by the integration interval.

For some applications, this accuracy is sufficient; for others, a better approximation is needed. The Euler method algorithm follows these steps:

1. Start with $t = t_{0}$, $x = x_{0}$.
2. Discretize time into equally spaced steps with spacing $h$, denoting each time point as $t_{i}$.
3. For each time point, find $x$ using the result of the previous iteration: $x_{i} = x_{i-1} + hf(x_{i-1}, t_{i-1})$.

## Runge-Kutta Method of Order 2

The RK2 method improves upon the Euler method by using the midpoint for evaluation, leading to better accuracy for the same step size $h$. Instead of evaluating the derivative at point $t$, as in Euler's method, RK2 evaluates it at the midpoint $t + h/2$.

By applying the Taylor series around the midpoint $t + h/2$, we derive:

$$
x(t + h) = x\left(t + \frac{h}{2}\right) + \frac{h}{2}\left(\frac{{\rm d}x}{{\rm d}t}\right)_{t+h/2} + \frac{h^2}{8}\left(\frac{{\rm d}^2x}{{\rm d}t^2}\right)_{t+h/2} + \mathcal{O}(h^3).
$$

Similarly, for $x(t)$, we have:

$$
x(t) = x\left(t + \frac{h}{2}\right) - \frac{h}{2}\left(\frac{dx}{dt}\right)_{t+h/2} + \frac{h^2}{8}\left(\frac{d^2x}{dt^2}\right)_{t+h/2} + \mathcal{O}(h^3).
$$

Subtracting these equations gives:

$$
x(t + h) = x(t) + h\left(\frac{dx}{dt}\right)_{t+h/2} + \mathcal{O}(h^3).
$$

Thus, the RK2 method is:

$$
\boxed{x(t + h) = x(t) + hf\left(x(t + h/2), t + h/2\right) + \mathcal{O}(h^3)}.
$$

The second-order term disappears, and the error is $\mathcal{O}(h^3)$, enhancing computational efficiency.

However, we need the function value at the midpoint $x(t + h/2)$, which is approximated using Euler's method with a half-step $h/2$:

$$
x(t + h/2) \approx x(t) + \frac{h}{2}f(x,t).
$$

Therefore, RK2 equations are:

1. Compute $k_{1} = hf(x,t)$.
2. Compute $k_{2} = hf\left(x + \frac{k_{1}}{2}, t + \frac{h}{2}\right)$.
3. Update $x(t + h) = x(t) + k_{2}$.

The local approximation error per step is $\mathcal{O}(h^3)$, while the global error is $\mathcal{O}(h^2)$, ensuring a more accurate solution than the Euler method for the same step size.

## Runge-Kutta Method of Order 4

The RK4 method extends RK2 by using more points between $x(t)$ and $x(t + h)$ through Taylor expansions, providing a better balance between complexity and accuracy. It is the most widely used method for solving ODEs due to its efficiency and precision. For most applications, the RK4 method is the de facto standard due to its ease of implementation and accurate results. The local approximation error per step is $\mathcal{O}(h^5)$, while the global error is approximately $\mathcal{O}(h^4)$, providing a precise and efficient solution.

We can set the RK4 algorithm such that:

1. **Calculate the slopes:**
    - $k_{1} = hf(x, t)$
    - $k_{2} = hf\left(x + \frac{k_{1}}{2}, t + \frac{h}{2}\right)$
    - $k_{3} = hf\left(x + \frac{k_{2}}{2}, t + \frac{h}{2}\right)$
    - $k_{4} = hf\left(x + k_{3}, t + h\right)$

2. **Update the function value:**
    - $x(t + h) = x(t) + \frac{1}{6}(k_{1} + 2k_{2} + 2k_{3} + k_{4})$

This method balances computational efficiency and accuracy, making it a preferred choice for solving differential equations in various scientific and engineering applications.

## Project layout

    
    mkdocs.yml    # The configuration file.
    docs/
	index.md  # The documentation homepage.
    ODE/
        ODE.py   # The code documented 
        __init__.py # Little description

