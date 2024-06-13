### This project includes two main questions about neural network controllers and fuzzy PID controllers
## Question 1: Neural Network Controller

  # Neural Network Structure as Controller:
  The task involves using a neural network directly as a controller for a nonlinear time-varying plant.
  It details the update rules for the neural network controller parameters to minimize an objective function using the chain rule.
  It includes mathematical formulations for the forward phase, calculation of control signals u(k)u(k), and update rules for the weights in both the output and input layers.
  The final weights and biases are updated iteratively to minimize the cost function.

  # MATLAB Implementation:
  The assignment requires writing MATLAB code to implement an adaptive neural network controller using gradient descent to minimize the objective function.
  The initial values of the weights are provided, and the reference signal is given in Figure 2.
  The results should include plots showing the alternation of controller parameters, reference-system signal, and control signal over 100 iterations.

# Question 2: PID Type Fuzzy Controller

  # Fuzzy PID Controller Design:
  The problem involves designing a PID-type fuzzy controller using triangular membership functions.
  Steps for the control procedure include tracking error calculation, obtaining controller inputs, determining the fuzzy controller's output, calculating the current control signal, and applying it to the system.

  # MATLAB Implementation:
  The assignment includes designing a fuzzy PID controller in MATLAB to minimize the tracking error for a given reference signal.

  # Proposing an Adaptation Mechanism:
  The final part involves proposing an adaptation mechanism for the controller based on neural networks, fuzzy logic, or another intelligent method.
  The Adaptive Neuro-Fuzzy Inference System (ANFIS) is a potential method, that integrates neural networks and fuzzy logic principles.
  ANFIS architecture is described, including the premise and consequence parts, and the five layers: fuzzification, rule, normalization, defuzzification, and output.
