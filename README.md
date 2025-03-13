# MST-CFD
一个“适应非结构网格的（MSH格式），基于有限体积的，求解单相单组分的可压层流问题的，二维三维通用的”CFD程序包，包含MSH文件格式读取程序，包含二阶TVD密度基求解器，面通量通过ROE或AUSM格式获取，粑粑山。
添加了非结构网格下可压瞬态SIMPLE、COUPLED算法，但是计算效果不好，计算SOD管结果还可以，但是任何更复杂模型都会发散。

A "unstructured mesh-adapted (for MSH), finite volume method - based, two-dimensional & three-dimensional general purpose CFD package for solving single-phase single-component laminary-flow problems", with .MSH file  reader, with second-order TVD density-based solver, with surface flux obtained by ROE or AUSM format, Pop-mount.
The compressible transient SIMPLE and COUPLED algorithms under unstructured grids have been added, but the calculation results are not satisfactory. The results of the SOD tube calculation are acceptable, but any other computation with more complex model will diverge.


