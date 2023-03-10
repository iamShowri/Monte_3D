function [P_mat, R_vec] = m3d_datastruct(n_total, box_length, R_gen, sortstr)

  %% Allocate memory:
  P_mat(n_total,4) = 0;
  R_vec(n_total,2) = 0;

  %% Assign particles to vertices of domain:
  P_mat(1,:) = [1, box_length, box_length, box_length];
  P_mat(2,:) = [2, 2*box_length, box_length, box_length];
  P_mat(3,:) = [3, 2*box_length, 2*box_length, box_length];
  P_mat(4,:) = [4, box_length, 2*box_length, box_length];
  P_mat(5,:) = [5, box_length, box_length, 2*box_length];
  P_mat(6,:) = [6, 2*box_length, box_length, 2*box_length];
  P_mat(7,:) = [7, 2*box_length, 2*box_length, 2*box_length];
  P_mat(8,:) = [8, box_length, 2*box_length, 2*box_length];

  R_vec(1:8,1) = (1:8);
  if (strcmp(sortstr, 'n') == 1)
    R_vec(1:8,2) = R_gen(1,2);
  else
    R_vec(1:8,2) = mean(R_gen(:,2));
  end

  %% Map data structure:
  P_mat((n_total+1):(n_total+8),2:4) = P_mat(1:8,2:4)-box_length;
  R_vec((n_total+1):(n_total+8),2) = R_vec(1:8,2);
  P_mat((2*n_total+1):(2*n_total+8),2) = P_mat(1:8,2);
  P_mat((2*n_total+1):(2*n_total+8),3:4) = P_mat(1:8,3:4)-box_length;
  R_vec((2*n_total+1):(2*n_total+8),2) = R_vec(1:8,2);
  P_mat((3*n_total+1):(3*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((3*n_total+1):(3*n_total+8),3:4) = P_mat(1:8,3:4)-box_length;
  R_vec((3*n_total+1):(3*n_total+8),2) = R_vec(1:8,2);
  P_mat((4*n_total+1):(4*n_total+8),2) = P_mat(1:8,2)-box_length;
  P_mat((4*n_total+1):(4*n_total+8),3) = P_mat(1:8,3);
  P_mat((4*n_total+1):(4*n_total+8),4) = P_mat(1:8,4)-box_length;
  R_vec((4*n_total+1):(4*n_total+8),2) = R_vec(1:8,2);
  P_mat((5*n_total+1):(5*n_total+8),2:3) = P_mat(1:8,2:3);
  P_mat((5*n_total+1):(5*n_total+8),4) = P_mat(1:8,4)-box_length;
  R_vec((5*n_total+1):(5*n_total+8),2) = R_vec(1:8,2);
  P_mat((6*n_total+1):(6*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((6*n_total+1):(6*n_total+8),3) = P_mat(1:8,3);
  P_mat((6*n_total+1):(6*n_total+8),4) = P_mat(1:8,4)-box_length;
  R_vec((6*n_total+1):(6*n_total+8),2) = R_vec(1:8,2);
  P_mat((7*n_total+1):(7*n_total+8),2) = P_mat(1:8,2)-box_length;
  P_mat((7*n_total+1):(7*n_total+8),3) = P_mat(1:8,3)+box_length;
  P_mat((7*n_total+1):(7*n_total+8),4) = P_mat(1:8,4)-box_length;
  R_vec((7*n_total+1):(7*n_total+8),2) = R_vec(1:8,2);
  P_mat((8*n_total+1):(8*n_total+8),2) = P_mat(1:8,2);
  P_mat((8*n_total+1):(8*n_total+8),3) = P_mat(1:8,3)+box_length;
  P_mat((8*n_total+1):(8*n_total+8),4) = P_mat(1:8,4)-box_length;
  R_vec((8*n_total+1):(8*n_total+8),2) = R_vec(1:8,2);
  P_mat((9*n_total+1):(9*n_total+8),2:3) = P_mat(1:8,2:3)+box_length;
  P_mat((9*n_total+1):(9*n_total+8),4) = P_mat(1:8,4)-box_length;
  R_vec((9*n_total+1):(9*n_total+8),2) = R_vec(1:8,2);
  P_mat((10*n_total+1):(10*n_total+8),2:3) = P_mat(1:8,2:3)-box_length;
  P_mat((10*n_total+1):(10*n_total+8),4) = P_mat(1:8,4);
  R_vec((10*n_total+1):(10*n_total+8),2) = R_vec(1:8,2);
  P_mat((11*n_total+1):(11*n_total+8),2) = P_mat(1:8,2);
  P_mat((11*n_total+1):(11*n_total+8),3) = P_mat(1:8,3)-box_length;
  P_mat((11*n_total+1):(11*n_total+8),4) = P_mat(1:8,4);
  R_vec((11*n_total+1):(11*n_total+8),2) = R_vec(1:8,2);
  P_mat((12*n_total+1):(12*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((12*n_total+1):(12*n_total+8),3) = P_mat(1:8,3)-box_length;
  P_mat((12*n_total+1):(12*n_total+8),4) = P_mat(1:8,4);
  R_vec((12*n_total+1):(12*n_total+8),2) = R_vec(1:8,2);
  P_mat((13*n_total+1):(13*n_total+8),2) = P_mat(1:8,2)-box_length;
  P_mat((13*n_total+1):(13*n_total+8),3:4) = P_mat(1:8,3:4);
  R_vec((13*n_total+1):(13*n_total+8),2) = R_vec(1:8,2);
  P_mat((14*n_total+1):(14*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((14*n_total+1):(14*n_total+8),3:4) = P_mat(1:8,3:4);
  R_vec((14*n_total+1):(14*n_total+8),2) = R_vec(1:8,2);
  P_mat((15*n_total+1):(15*n_total+8),2) = P_mat(1:8,2)-box_length;
  P_mat((15*n_total+1):(15*n_total+8),3) = P_mat(1:8,3)+box_length;
  P_mat((15*n_total+1):(15*n_total+8),4) = P_mat(1:8,4);
  R_vec((15*n_total+1):(15*n_total+8),2) = R_vec(1:8,2);
  P_mat((16*n_total+1):(16*n_total+8),2) = P_mat(1:8,2);
  P_mat((16*n_total+1):(16*n_total+8),3) = P_mat(1:8,3)+box_length;
  P_mat((16*n_total+1):(16*n_total+8),4) = P_mat(1:8,4);
  R_vec((16*n_total+1):(16*n_total+8),2) = R_vec(1:8,2);
  P_mat((17*n_total+1):(17*n_total+8),2:3) = P_mat(1:8,2:3)+box_length;
  P_mat((17*n_total+1):(17*n_total+8),4) = P_mat(1:8,4);
  R_vec((17*n_total+1):(17*n_total+8),2) = R_vec(1:8,2);
  P_mat((18*n_total+1):(18*n_total+8),2:3) = P_mat(1:8,2:3)-box_length;
  P_mat((18*n_total+1):(18*n_total+8),4) = P_mat(1:8,4)+box_length;
  R_vec((18*n_total+1):(18*n_total+8),2) = R_vec(1:8,2);
  P_mat((19*n_total+1):(19*n_total+8),2) = P_mat(1:8,2);
  P_mat((19*n_total+1):(19*n_total+8),3) = P_mat(1:8,3)-box_length;
  P_mat((19*n_total+1):(19*n_total+8),4) = P_mat(1:8,4)+box_length;
  R_vec((19*n_total+1):(19*n_total+8),2) = R_vec(1:8,2);
  P_mat((20*n_total+1):(20*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((20*n_total+1):(20*n_total+8),3) = P_mat(1:8,3)-box_length;
  P_mat((20*n_total+1):(20*n_total+8),4) = P_mat(1:8,4)+box_length;
  R_vec((20*n_total+1):(20*n_total+8),2) = R_vec(1:8,2);
  P_mat((21*n_total+1):(21*n_total+8),2) = P_mat(1:8,2)-box_length;
  P_mat((21*n_total+1):(21*n_total+8),3) = P_mat(1:8,3);
  P_mat((21*n_total+1):(21*n_total+8),4) = P_mat(1:8,4)+box_length;
  R_vec((21*n_total+1):(21*n_total+8),2) = R_vec(1:8,2);
  P_mat((22*n_total+1):(22*n_total+8),2:3) = P_mat(1:8,2:3);
  P_mat((22*n_total+1):(22*n_total+8),4) = P_mat(1:8,4)+box_length;
  R_vec((22*n_total+1):(22*n_total+8),2) = R_vec(1:8,2);
  P_mat((23*n_total+1):(23*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((23*n_total+1):(23*n_total+8),3) = P_mat(1:8,3);
  P_mat((23*n_total+1):(23*n_total+8),4) = P_mat(1:8,4)+box_length;
  R_vec((23*n_total+1):(23*n_total+8),2) = R_vec(1:8,2);
  P_mat((24*n_total+1):(24*n_total+8),2) = P_mat(1:8,2)-box_length;
  P_mat((24*n_total+1):(24*n_total+8),3:4) = P_mat(1:8,3:4)+box_length;
  R_vec((24*n_total+1):(24*n_total+8),2) = R_vec(1:8,2);
  P_mat((25*n_total+1):(25*n_total+8),2) = P_mat(1:8,2);
  P_mat((25*n_total+1):(25*n_total+8),3:4) = P_mat(1:8,3:4)+box_length;
  R_vec((25*n_total+1):(25*n_total+8),2) = R_vec(1:8,2);
  P_mat((26*n_total+1):(26*n_total+8),2) = P_mat(1:8,2)+box_length;
  P_mat((26*n_total+1):(26*n_total+8),3:4) = P_mat(1:8,3:4)+box_length;
  R_vec((26*n_total+1):(26*n_total+8),2) = R_vec(1:8,2);

end
