# Tools-for-processing-disulfide-bond-a-pymol-plugin
I used Gemini 3 to write a plugin for handling the disulfide bond deficiency issue after Swiss-Model homology modeling.This plugin is written in Python and relies on the PyMOL software. For more usage instructions, please refer to the README file.
SS-Manager: PyMOL Disulfide Bond Repair Suite
SS-Manager 是一个专为结构生物学和分子动力学（MD）模拟准备的全能 PyMOL 插件。

它旨在解决**同源建模（Homology Modeling）**中常见的痛点：

二硫键拓扑丢失： PDB 文件中只有原子坐标，没有 SSBOND 或 CONECT 记录。

严重的几何畸变： 建模补全的缺失环（Missing Loops）飘在溶剂中，导致二硫键距离过远（> 10Å），无法直接进行能量最小化。

格式兼容性差： 导出的 PDB 缺少 GROMACS/Amber 所需的 SSBOND 头信息。

📦 功能特性 (Features)
拓扑优先 (Topology First): 使用 CSV 文件作为二硫键的“真理来源 (Source of Truth)”，防止误连或漏连。

精确修复 (Surgical Repair): 提供 ss_snap 命令，通过刚体平移（Rigid Body Translation）精确修复断裂的大跨度二硫键。

诊断工具 (Diagnostics): 自动检测所有二硫键的键长，并标记状态（✅ OK / ⚠️ Stretched / 🚨 Broken）。

MD 友好 (Simulation Ready): 导出的 PDB 文件自动封装标准的 SSBOND 头信息，直接兼容 pdb2gmx 或 tleap。

拓扑克隆 (Topology Transfer): 将晶体结构（Template）的连接关系一键克隆给模型（Model）。

📥 安装 (Installation)
下载 ss_manager.py 到您的工作目录。

打开 PyMOL。

在 PyMOL 命令行或脚本中加载：

代码段

run ss_manager.py
🧪 快速上手 (Quick Start / Tutorial)
场景： 你有一个完美的晶体结构 (apo.pdb) 和一个同源建模得到的激活态结构 (activate.pdb)。activate 中有一个环是补全的，二硫键断开了 18Å，你需要修复它并准备跑 MD。

第一步：建立标准 (Define Topology)
从可靠的参考结构导出二硫键列表。

代码段

load apo.pdb
ss_export apo, bonds.csv
提示：打开生成的 bonds.csv，根据 UniProt 或文献手动检查，补全任何遗漏的键（格式：Chain1,Resi1,Chain2,Resi2,Note）。

第二步：诊断问题 (Diagnose)
加载目标模型，检查它的二硫键情况。

代码段

load activate.pdb
ss_check_dist activate, bonds.csv
输出示例：

A:670 - C:670 | 18.08 Å | BROKEN (Red) <- 这就是我们要修的！ A:22 - A:3 | 2.02 Å | OK (Green)

第三步：修复大裂谷 (Fix the Gap)
对于 > 5Å 的红色断裂键，使用 ss_snap 进行刚体对齐。 假设我们需要把 C 链的环（C:660-680）移过去找 A 链。

执行瞬移 (Snap):

代码段

# 语法: ss_snap [锚点原子], [移动原子], [移动的整体范围]
ss_snap /activate//A/670/SG, /activate//C/670/SG, /activate//C/660-680
缝合骨架 (Heal Backbone): 瞬移会导致骨架断裂，必须用 PyMOL 的 Sculpting 功能修复。

代码段

# 1. 保护不需要动的部分 (Fix Core)
select core, not (chain C and resi 660-680)
flag fix, core

# 2. 开启力场自动缝合
sculpt_activate all
# ... 等待 10 秒，看着断口自动愈合 ...
sculpt_deactivate all
flag free, all
第四步：建立连接并导出 (Finalize)
代码段

# 1. 在 PyMOL 中建立连线 (Visual)
ss_import activate, bonds.csv

# 2. 导出最终 PDB (自动写入 SSBOND)
ss_save_pdb activate, bonds.csv, activate_final.pdb
现在，activate_final.pdb 已经准备好进入 GROMACS 进行 pdb2gmx 了！

📚 命令参考手册 (Command Reference)
1. ss_export
导出二硫键列表到 CSV 文件。

用法: ss_export [object], [filename]

示例: ss_export 7xgd, my_bonds.csv

说明: 默认扫描距离 < 3.2Å 的 CYS-CYS 对。

2. ss_import
读取 CSV 并强制在 PyMOL 中建立 Bond 连接。

用法: ss_import [object], [filename]

说明: 这只是建立拓扑（画线），不会移动原子位置。距离过远的键会显示为粉色。

3. ss_check_dist
[推荐] 读取 CSV 并测量目标结构中的实际距离。

用法: ss_check_dist [object], [filename]

状态码:

🟢 OK (< 2.5Å): 完美。

🟡 STRETCHED (2.5 - 4.5Å): 拉伸。通常不需要手动修，MD 能量最小化可自动解决。

🔴 BROKEN (> 4.5Å): 断裂。必须使用 ss_snap 修复。

4. ss_snap
[核心功能] 执行精确的向量平移（Vector Translation），修复大跨度断裂。

用法: ss_snap [fixed_atom], [moving_atom], [moving_scope]

参数:

fixed_atom: 不动的锚点（目标位置）。

moving_atom: 需要移动的那个原子（当前位置）。

moving_scope: 需要跟随平移的整个原子选区（通常是一个 Loop）。

注意: 使用后必须配合 PyMOL 的 sculpt_activate 来消除原子碰撞。

5. ss_save_pdb
保存 PDB 文件并注入 SSBOND 记录。

用法: ss_save_pdb [object], [csv_source], [output_filename]

说明: 它是 PyMOL save 命令的封装，解决了 PyMOL 不写文件头的 Bug。

6. ss_compare (Legacy V11)
比较参考结构（Ref）与目标结构（Target）的二硫键差异。

用法: ss_compare [ref_obj], [target_obj]

说明: 使用序列相似度（Sequence Similarity）算法进行配对，能区分同源建模中的保守键和变异键。

7. ss_transfer (Legacy V12)
将 Source 的拓扑直接克隆给 Target。

用法: ss_transfer [source], [target]

说明: 适用于两个结构序列完全一致，且不需要 CSV 中间步骤的场景。

⚠️ 常见问题 (FAQ)
Q: 为什么 sculpt_activate 没反应？ A: 可能是距离太远（>10Å），超出了 PyMOL 简易力场的截断值。

解决: 先用 ss_snap 把它们瞬间拉近，然后再用 Sculpting。

Q: 粉色 (Magenta) 的键需要修吗？ A: 如果距离在 3.0 - 4.5 Å 之间，不需要。这属于轻微拉伸，GROMACS/Amber 的能量最小化（EM）步骤会根据力场参数自动把它拉回 2.05 Å。只有红色的（> 5 Å）才建议手动修。

Q: 为什么生成的 PDB 还是报错 "Missing atoms"？ A: 确保你在 ss_snap 后运行了 sculpting 来修复骨架断裂。刚体平移会撕裂肽链（C-N 键断开），Sculpting 会重新连接它们。

📝 License
MIT License. Feel free to use and modify for your research.
