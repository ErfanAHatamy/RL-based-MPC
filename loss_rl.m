function loss = loss_rl(s, a)

error = [s-[1;0]; a];

loss = error' * diag([1, 0.01, 0.01]) * error;

end


