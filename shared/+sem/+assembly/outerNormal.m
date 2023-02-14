function outerNormal(be, vxy, etov)
    let n = be.n
    let k = try! be.toEdgeNumber(etov: etov)
    
    let (n1,n2,n3) = etov[n]
    
    var (x1,y1): (Double, Double)
    var (x2,y2): (Double, Double)
    
    switch k {
    case .edge1:
        (x1,y1) = vxy[n1]
        (x2,y2) = vxy[n2]
    case .edge2:
        (x1,y1) = vxy[n2]
        (x2,y2) = vxy[n3]
    case .edge3:
        (x1,y1) = vxy[n3]
        (x2,y2) = vxy[n1]
    }
    
    let DeltaX = x2-x1
    let DeltaY = y2-y1
    
    let normalX = DeltaY/(sqrt(DeltaX*DeltaX + DeltaY*DeltaY))
    let normalY = -DeltaX/(sqrt(DeltaX*DeltaX + DeltaY*DeltaY))
    
    return (normalX,normalY)
}