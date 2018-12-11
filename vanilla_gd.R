library(pracma)
fun = function(x){
    return(1/2*x[1]^2+1/2*sin(3*x)+x[2]^2)
}

dfun = function(x){
    c(x[1]+3/2*cos(3*x), 2*x[2])
}

train_2d = function(trainer, params){
    results = c(params[1], params[2])
    for(i in 1:20){
        r = trainer(params)
        params[1] = r[1]; params[2] = r[2]
        message = paste("epoch", i, ",", params[1], ",", params[2])
        results = rbind(results, c(params[1], params[2]))
    }
    return(results)
}

showTrace = function(results){
    range_x = 3
    range_y = 3
    x <- seq(-range_x, range_x,length.out=100)
    y <- seq(-range_y, range_y,length.out=100)
    z <- outer(1/2*x^2 + 1/2*sin(3*x), y^2, `+`)
    view = contour(x, y, z, col = 'blue', nlevels = 100, drawlabel = FALSE, main='test')
    points(results[,1], results[,2], pch=20)
    lines(results[,1], results[,2])
}
# 1/2*x^2 + 1/2*sin(3*x)

vanilla_gd = function(params){
    x1 = params[1]
    y1 = params[2]
    eta = params[3]
    return(c(  x1-eta*(x1+3/2*cos(3*x1)), y1-eta*(2*y1)   ))
}

momentum = function(params){
    x1 = params[1]
    y1 = params[2]
    v1 = 0; v2 = 0;
    eta = params[3]
    gamma = params[4]
    v1 = gamma*v1 + eta*(x1+3/2*cos(3*x1))
    v2 = gamma*v2 + eta*(2*y1)
    return(c(x1-v1, y1-v2, v1, v2))
}

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



#persp(x,y,z, col = 'blue')

# (result = train_2d(vanilla_gd, c(0, 3, 0.1)))
#(result = train_2d(momentum, c(0, 3, 0.3, 0.5)))
#(result = train_2d(Newton, c(-1, 1)))
# (result = train_2d(dampedNewton, c(-1, 0.5)))


showTrace(result)







