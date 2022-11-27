function loss = loss_rl(s, s_r, a, a_r)

error = [s-s_r; a-a_r];

loss = error' * diag([1, 0.01, 0.01]) * error;

end


