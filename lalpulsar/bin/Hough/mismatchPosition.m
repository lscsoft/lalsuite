function retval = mismatchPosition(mismatch, lambda, alpha, delta)

  %% Rotation matrix
  M(1,1) = sin(alpha);
  M(1,2) = sin(delta)*cos(alpha);
  M(1,3) = cos(delta)*cos(alpha);
  M(2,1) = -cos(alpha);
  M(2,2) = sin(delta)*sin(alpha);
  M(2,3) = cos(delta)*sin(alpha);
  M(3,1) = 0;
  M(3,2) = -cos(delta);
  M(3,3) = sin(delta);

  %% input vector
  x(1) = sin(mismatch)*cos(lambda);
  x(2) = sin(mismatch)*sin(lambda);
  x(3) = cos(lambda);


  %% output vector
  y = M*x';

  %% convert y to spherical polar coordinates
  retval(1) = atan(y(2)/y(1));
  retval(2) = asin(y(3));

endfunction

