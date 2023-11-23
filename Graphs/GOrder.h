//
// Created by ldd on 2023/6/18.
//

#ifndef MLCDEC_GORDER_H
#define MLCDEC_GORDER_H


enum G_ORDER {
    DEFAULT, DEN_INC, DEN_DEC, CORE_INC, CORE_DEC, E_INC, E_DEC
};

static string GraphOrder2Str(G_ORDER order) {
    if (order == DEFAULT) return "rd";
    else if (order == DEN_INC) return "di";
    else if (order == DEN_DEC) return "dd";
    else if (order == CORE_INC) return "ci";
    else if (order == CORE_DEC) return "cd";
    else if (order == E_INC) return "ei";
    else if (order == E_DEC) return "ed";
    return "";
}

#endif //MLCDEC_GORDER_H


