for f = 1:n
    
    j1 = face(f,2);
    j2 = face(f,3);
    j3 = face(f,4);
    x = [vertex(j1,2);vertex(j2,2);vertex(j3,2)];
    y = [vertex(j1,3);vertex(j2,3);vertex(j3,3)];
    z = [vertex(j1,4);vertex(j2,4);vertex(j3,4)];
    
    h = patch(x,y,z,'c');
    set(h,'edgecolor','r')

end