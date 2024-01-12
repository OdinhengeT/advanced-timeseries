function y_t1 = EEstep(f, param, t0, y_t0, dt)  
    y_t1 = y_t0 + dt .* f(t0, y_t0, param);
end

