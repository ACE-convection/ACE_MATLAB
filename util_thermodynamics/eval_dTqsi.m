function dTqsi = eval_dTqsi(T,p,qt)
%{
  Calculate partial derivative of (ice) saturation specific humidity with respect to temperature.

  Units: kg/kg/T.

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    qt (double): total water mixing ratio in kg/kg.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.

  Returns:
    double: derivative of saturation specific humidity w.r.t. T.
%}
  esi = eval_esi(T);
  % With es from Eq. (4.4.15) in Emanual (1994)
  %   accurate within <1e-3 rel. error for -30<=T<=0-deg-C
  % Rel. error increases exponentially with T
  if nargin==2 %assuming condensate-free
    dTqsi = 0.621971830985916.*p./( p-esi*(1-0.621971830985916) ).^2 .* (6111.72784./T+0.15215).*esi./T;
	else %given qt (suitable for plume computations)
    dTqsi = 0.621971830985916.*p.*(1-qt)./(p-esi).^2 .* (6111.72784./T+0.15215).*esi./T;
  end