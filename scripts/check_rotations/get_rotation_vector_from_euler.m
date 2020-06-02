function v = get_rotation_vector_from_euler( phi, theta, psi );
v = SpinCalc( 'EA313toEV', [phi,theta,psi],1.0e-5,0 );
v(:,4) = v(:,4) - 360*(v(:,4)>180);
v = v(:,[1:3]) .* repmat( v(:,4),[1 3] );   
    
