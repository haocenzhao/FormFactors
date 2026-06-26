# print_gegm_table.py
# 调用已有的 GetFF(kID, kQ2)，输出表格中各 Q2 对应的中心值

from GetFF import GetFF

# 左侧表：mu_p * GEp / GMp
Q2_p_list = [3.4, 5.5, 7.8, 11.7]

# 右侧表：mu_n * GEn / GMn
Q2_n_list = [2.9, 6.6, 9.9, 4.5]


def calc_mu_p_ge_over_gm(Q2):
    gep_red, _ = GetFF(1, Q2)   # GEp / (GEp(0)*GD), 其中 GEp(0)=1
    gmp_red, _ = GetFF(2, Q2)   # GMp / (mu_p*GD)
    return gep_red / gmp_red    # = mu_p * GEp / GMp


def calc_mu_n_ge_over_gm(Q2):
    gen_red, _ = GetFF(3, Q2)   # GEn / (GEn(0)*GD)  按你当前代码约定
    gmn_red, _ = GetFF(4, Q2)   # GMn / (mu_n*GD)
    return gen_red / gmn_red    # = mu_n * GEn / GMn


def main():
    print("=" * 78)
    print(f"{'Q2(p)':>8}  {'mu_p*GEp/GMp':>18}    {'Q2(n)':>8}  {'mu_n*GEn/GMn':>18}")
    print("=" * 78)

    nrows = max(len(Q2_p_list), len(Q2_n_list))

    for i in range(nrows):
        if i < len(Q2_p_list):
            q2p = Q2_p_list[i]
            val_p = calc_mu_p_ge_over_gm(q2p)
            left = f"{q2p:8.1f}  {val_p:18.8f}"
        else:
            left = f"{'':>8}  {'':>18}"

        if i < len(Q2_n_list):
            q2n = Q2_n_list[i]
            val_n = calc_mu_n_ge_over_gm(q2n)
            right = f"{q2n:8.1f}  {val_n:18.8f}"
        else:
            right = f"{'':>8}  {'':>18}"

        print(left + "    " + right)

    print("=" * 78)


if __name__ == "__main__":
    main()