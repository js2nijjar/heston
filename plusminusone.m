function result = plusminusone()
    r = rand;
    threshold = 0.5;
    if r < threshold
        result = -1;
    else
        result = 1;
    end
end