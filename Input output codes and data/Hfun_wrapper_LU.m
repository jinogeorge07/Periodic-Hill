function y = Hfun_wrapper_LU(x, tflag, Aop, Aadjop, C)
%Hfun_wrapper_LU   Operator for ARPACK with cached LU factors
%
%  y = Hfun_wrapper_LU(x, 'notransp', Aop, Aadjop, C)   returns C*(M\ (Btilde*x))
%  y = Hfun_wrapper_LU(x, 'transp'  , Aop, Aadjop, C)   returns    Aadjop(x)

  switch tflag
    case 'notransp'
      % forward action: y = C * (M \ (B_tilde * x))
      y = C * Aop(x);

    case 'transp'
      % adjoint   : y = M' \ (C' * x)  then times Btilde'
      % but since Aadjop = M' \ (C' * x) *and* Btilde' is built into that,
      % we just call:
      y = Aadjop(x);

    otherwise
      error('Hfun_wrapper_LU: unknown tflag %s', tflag);
  end
end