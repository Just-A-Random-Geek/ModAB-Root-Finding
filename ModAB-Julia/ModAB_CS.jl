function modab_CS(f, left::Real, right::Real, target::Real=0.0; precision::Float64=1e-14, maxIter::Int=200)
    x1, x2 = min(left, right), max(left, right)
    y1 = f(x1) - target
    abs(y1) <= precision && return x1
    y2 = f(x2) - target
    abs(y2) <= precision && return x2
    eps1 = precision * 1e-3
    eps2 = precision * (x2 - x1)
    if abs(target) >= 1
        eps1 *= abs(target)
    else
        eps1 = 0
    end
    side = 0
    bisection = true
    C = 16 # safety factor for threshold corresponding to 4 iterations = 2^4
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    # calculate k on each bisection step with account for local function properties and symmetry
    for i in 1:maxIter
        local x3, y3
        if bisection
            x3 = (x1 + x2) / 2
            y3 = f(x3) - target  # Function value at midpoint
            ym = (y1 + y2) / 2 # Ordinate of chord at midpoint
            r = 1 - abs(ym / (y2 - y1)) # Symmetry factor
            k = r * r # Deviation factor
            # Check if the function is close enough to linear
            if abs(ym - y3) < k * (abs(y3) + abs(ym))
                bisection = false
                threshold = (x2 - x1) * C
            end
        else
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
            if x3 <= x1
                x3 = x1
                y3 = y1
            elseif x3 >= x2
                x3 = x2
                y3 = y2
            else
                y3 = f(x3) - target
            end    
            threshold /= 2
        end

        if abs(y3) <= eps1 || abs(x2 - x1) <= eps2 # Convergence check
            return x3
        end

        if sign(y1) == sign(y3)
            if side == 1
                m = 1 - y3 / y1
                if m <= 0
                    y2 /= 2
                else
                    y2 *= m
                end
            elseif !bisection
                side = 1
            end
            x1, y1 = x3, y3
        else
            if side == -1
                m = 1 - y3 / y2
                if m <= 0
                    y1 /= 2
                else
                    y1 *= m
                end
            elseif !bisection
                side = -1
            end
            x2, y2 = x3, y3
        end
        if x2 - x1 > threshold # AB failed to shrink the interval enough
            bisection = true
            side = 0
        end
    end
    return NaN
end