library(pracma) # for line search
fun = function(x){
    return(1/2*x[1]^2+1/2*sin(3*x[1])+x[2]^2)
}

dfun = function(x){
    return(c(x[1]+3/2*cos(3*x[1]), 2*x[2]))
}

train_2d = function(trainer, params, epoch=20){
    results = c(params[1], params[2])
    for(i in 1:epoch){
        params = trainer(params)
        if(params[1] == 0) break
        message = paste("epoch", i, ",", params[1], ",", params[2])
        results = rbind(results, c(params[1], params[2]))
    }
    return(results)
}

showTrace = function(results, title_){
    range_x = 2
    range_y = 2
    x <- seq(-range_x, range_x,length.out=100)
    y <- seq(-range_y, range_y,length.out=100)
    z <- outer(1/2*x^2 + 1/2*sin(3*x), y^2, `+`)
    view = contour(x, y, z, col = 'blue', nlevels = 40, drawlabel = FALSE, main=title_)
    points(results[,1], results[,2], pch=20)
    lines(results[,1], results[,2])
}

addTrace = function(result, color){
    points(result[,1], result[,2], pch=20,col=color)
    lines(result[,1], result[,2], col=color)
}

objLineChart = function(result, type_, color){
    y = c()
    epoch = 1:nrow(result)
    for(i in epoch){
        y = c(y, fun(result[i, ]))
    }
    if(type_ == 1){
        plot(epoch, y, pch = 20, type="o", col=color)
    }else{
        lines(epoch, y, pch=20, type="o", col=color)   
    }
}

addLine = function(result){
    y = c()
    x = 1:nrow(result)
    for(i in x){
        y = c(y, fun(result[i, ]))
    }
}

# Optimization functions
vanilla_gd = function(params){
    x1 = params[1]
    y1 = params[2]
    x = c(params[1], params[2])
    eta = params[3]
    x = x - eta*dfun(x)
    return(c(x, eta))
}

momentum = function(params){
    x = c(params[1], params[2])
    v = c(params[5], params[6])
    eta = params[3]
    gamma = params[4]
    v = gamma*v + eta*dfun(x)
    x = x-v
    return(c(x, eta, gamma, v))
}

NAG = function(params){
    x = c(params[1], params[2])
    v = c(params[5], params[6])
    eta = params[3]
    gamma = params[4]
    v = gamma*v + eta*dfun(x - gamma*v)
    x = x-v
    return(c(x, eta, gamma, v))
}

# 基本牛顿法
Newton = function(params){
    x = params[1]
    y = params[2]
    gx = x + 3/2*cos(3*x)
    gy = 2*y
    gxx = 1 - 9/2*sin(3*x)
    gyy = 2
    return(c(x - gx/gxx, y - gy/gyy))
}

dampedNewton = function(params){
    x = params[1]
    y = params[2]
    gx = x + 3/2*cos(3*x)
    gy = 2*y
    gxx = 1 - 9/2*sin(3*x)
    gyy = 2
    lambda = softline(c(x,y), c(-gx/gxx, -gy/gyy), fun, dfun)
    return(c(x,y) + lambda*c(-gx/gxx, -gy/gyy))
}

adagrad = function(params){
    x = c(params[1], params[2])
    s = c(params[3], params[4])
    eta = params[5]
    eps = 1e-6
    g = dfun(x)
    s = s + g*g
    x = x - eta/sqrt(s+eps)*g
    return(c(x, s, eta))
}


rmsprop = function(params){
    x = c(params[1], params[2])
    s = c(params[3], params[4])
    eta = params[5]
    gamma = params[6]
    eps = 1e-8
    g = dfun(x)
    s = gamma*s + (1-gamma)*g*g
    x = x - eta/sqrt(s + eps)*g
    return(c(x, s, eta, gamma))
}

Adam = function(params){
    x = c(params[1], params[2])
    s = c(params[3], params[4])
    eta = params[5]
    v = c(params[6], params[7])
    beta1 = params[8]
    beta2 = params[9]
    t = params[10]
    eps = 1e-8
    g = dfun(x)
    v = beta1*v + (1-beta1)*g
    s = beta2*s + (1-beta2)*g*g
    v_hat = v / (1-beta1^t)
    # print(paste('s', s))
    s_hat = s / (1-beta2^t)
    # print(paste('s_hat', s_hat))
    g_t = eta * v_hat / sqrt(s_hat + eps)
    x_t = x - g_t
    t = t + 1
    # print(c(x_t, s_hat, eta, v_hat, beta1, beta2, t))
    return(c(x_t, s_hat, eta, v_hat, beta1, beta2, t))
}

conjugateGradient = function(x, eps){
    results = x
    r = 20;
    for(k in 0:r){
        if(k == 0){
            p = - dfun(x)
            g = dfun(x)
        }
        lambda = softline(x, p, fun, dfun)
        x = x + lambda*p
        g0 = g
        g = dfun(x)
        if(sqrt(g*g) < eps) break
        if(k==r-1){
            k = 0
        }else{
            alpha = g*g/(g0*g0)
            p = -g + alpha*p
        }
        if(p %*% g){
            k = 0
        }
        results = rbind(results, x)
    }
    return(results)
}

# Examples
# (result = train_2d(vanilla_gd, c(0.8, -1, 0.1)))
# (result = train_2d(momentum, c(0.5, 1, 0.3, 0.5, 0, 0)))
# (result = train_2d(NAG, c(0.5, 1.5, 0.1, 0.6, 0, 0)))
# (result = train_2d(Newton, c(-1, 0.5)))
# (result = train_2d(Newton, c(-0.8, 0.5)))
# (result = train_2d(dampedNewton, c(-1, 0.8)))
# (result = train_2d(adagrad, c(-1, 0.5, 0, 0, 0.1)))
# (result = train_2d(rmsprop, c(-1, 0.5, 0, 0, 0.1, 0.9)))
# (result = train_2d(Adam, c(-1, 0.5, 0, 0, 0.1, 0, 0, 0.9, 0.999, 1)))
# (result = conjugateGradient(c(-1, 1), 0.000001))

# showTrace(result, 'vanilla')






