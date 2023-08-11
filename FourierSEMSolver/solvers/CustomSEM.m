function customSEM = CustomSEM(Mx, Sx, V, F_filter, conn, etov, P, interface)
    customSEM = struct('Mx', Mx, 'Sx', Sx, 'V', V, ...\
        'F_filter', F_filter, 'conn', conn, 'etov', etov, ...\
        'scheme_order', P, 'interface', interface);
end