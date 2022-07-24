======================================================
``odeint``: Ordinary Differential Equations Integrator
======================================================

The integration of the equations of motions is taken care by the ``odeint`` module. It is based on `boost::numeric::odeint <https://www.boost.org/doc/libs/1_66_0/libs/numeric/odeint/doc/html/index.html>`_ but restricted to the needs of the project. 

How does it work?
-----------------

The problems are described in the form

.. math::
	\begin{align}
		x'(t) &= f(x, t)\\
		x(0) &= x_0.
	\end{align}

The numerical integration methods calculate the solution to :math:`x(t)` iteratively, in general there is a formula that obtains :math:`x(t + \Delta t)` as a function of :math:`x(t)`, :math:`f` and :math:`\Delta t`. For example, the Euler method approximates 

.. math::
	x(t + \Delta t) = x(t) + x'(t) \Delta t.

In ``odeint`` the methods of integration are implemented in *stepper* classes that are only required to have a method ``do_step(system, state, current_time, delta_time)`` that takes the equation ``system`` :math:`f`, the current ``state`` :math:`x(t)`, current time :math:`t` and the time delta :math:`\Delta t`, then applies the method in order to approximate the state :math:`x(t + \Delta t)` and outputs to the out parameter ``state``.

The ``do_step`` method is called iteratively in order to integrate the equation to the desired time using the following function:

.. doxygenfunction:: integrate(stepper_type& stepper, system_type& sys, state_type& x, scalar_type t0, scalar_type dt, size_t Nsteps, observer_type& obs, size_t obs_skip_steps = 0)

This function accepts templated arguments so that it can use different integration methods, equation systems, and observers. The function loops ``Nsteps`` times and evaluates the ``stepper.do_step()`` with ``sys``, the current time and state and then updates the state and the time (``+= dt``). ``obs`` is called every ``obs_skip_steps + 1`` to observe (output to file or save) the current state.

.. note::
	In the future it might be useful or necessary to implement other integration function that uses adaptative time steps.

Implemented steppers
--------------------

The currently implemented steppers are 

.. doxygenclass:: EulerStepper
.. doxygenclass:: RK46NL