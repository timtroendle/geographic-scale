localrules: figures


rule figures:
    message: "Collect and rename all figures."
    input:
        "build/output/{resolution}/cost.{plot_suffix}",
        "build/output/{resolution}/cost-uncertainty.{plot_suffix}",
        "build/output/{resolution}/composition.{plot_suffix}",
        "build/output/{resolution}/composition-uncertainty.{plot_suffix}",
        "build/output/{resolution}/map-energy.{plot_suffix}",
        "build/output/{resolution}/map-cost.{plot_suffix}",
        "build/output/{resolution}/total-sobol-diff.{plot_suffix}"
    output:
        "build/output/{resolution}/figures/figure1.{plot_suffix}",
        "build/output/{resolution}/figures/figure2.{plot_suffix}",
        "build/output/{resolution}/figures/figure3.{plot_suffix}",
        "build/output/{resolution}/figures/figure4.{plot_suffix}",
        "build/output/{resolution}/figures/figure5.{plot_suffix}",
        "build/output/{resolution}/figures/figure6.{plot_suffix}",
        "build/output/{resolution}/figures/figure7.{plot_suffix}"
    run:
        from shutil import copyfile
        for i in range(len(output)):
            copyfile(input[i], output[i])
