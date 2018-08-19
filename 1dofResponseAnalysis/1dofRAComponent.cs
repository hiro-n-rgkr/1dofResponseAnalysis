/**
Copyright (c) 2018 hiron_rgrk
This software is released under the MIT License.
See LICENSE
**/

using System;
using System.Collections.Generic;

using Grasshopper.Kernel;

using Rhino;
using Rhino.Geometry;


/// <summary>
/// コンポーネントの定義
/// </summary>
namespace GH_NewmarkBeta
{
    public class NewmarkBetaComponet : GH_Component
    {
        public NewmarkBetaComponet()
            : base("1dof Response Analysis",                // 名称
                   "1dof RA",                               // 略称
                   "Response Analysis of the Single dof",   // コンポーネントの説明
                   "rgkr",                                  // カテゴリ(タブの表示名)
                   "Response Analysis"                      // サブカテゴリ(タブ内の表示名)
                  )
        {
        }

        /// <summary>
        /// インプットパラメータの登録
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Mass", "M", "Lumped Mass(ton)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Stiffness", "K", "Spring Stiffness(kN/m)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Damping ratio", "h", "Damping ratio", GH_ParamAccess.item);
            pManager.AddNumberParameter("Time Increment", "dt", "Time Increment(sec)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Beta", "Beta", "Parameters of Newmark β ", GH_ParamAccess.item, 0.25);
            pManager.AddIntegerParameter("N", "N", "Parameters of Newmark β ", GH_ParamAccess.item);
            pManager.AddTextParameter("WAVE", "WAVE", "Parameters of Newmark β ", GH_ParamAccess.item);
        }

        /// <summary>
        /// アウトプットパラメータの登録
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Model", "Model", "output Model Data", GH_ParamAccess.item);
            pManager.AddNumberParameter("Acceleration", "Acc", "output Acceleration(cm/s^2)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Velocity", "Vel", "output Velocity(cm/s)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Displacement", "Disp", "output Displacement(cm)", GH_ParamAccess.item);
        }

        /// <summary>
        /// 解析を実行する箇所
        /// </summary>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // パラメータの定義 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            List<double> Model = new List<double>();
            double M = double.NaN;    // 質量
            double K = double.NaN;    // 剛性
            double h = double.NaN;    // 減衰定数
            double dt = double.NaN;   // 時間刻み
            double beta = double.NaN; // 解析パラメータ
            int N = 0;                // 波形データ数
            string wave_str = "0";

            // grasshopper からデータ取得　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            if (!DA.GetData(0, ref M)) { return; }
            if (!DA.GetData(1, ref K)) { return; }
            if (!DA.GetData(2, ref h)) { return; }
            if (!DA.GetData(3, ref dt)) { return; }
            if (!DA.GetData(4, ref beta)) { return; }
            if (!DA.GetData(5, ref N)) { return; }
            if (!DA.GetData(6, ref wave_str)) { return; }

            // モデルデータ出力用＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            Model.Add(M);
            Model.Add(K);

            //　地震波データの処理　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            //　カンマ区切りで波形を入力するので、カンマで区切り配列に入れている
            char[] delimiter = { ',' };    //分割文字
            double[] wave = new double[N];
            string[] wk;
            wk = wave_str.Split(delimiter);  //カンマで分割
            for (int i = 0; i < N; i++)
            {
                wave[i] = double.Parse(wk[i]);
            }

            //　応答解析　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            double[] out_a = new double[N];
            double[] out_v = new double[N];
            double[] out_d = new double[N];

            Solver.NewmarkBeta_solver slv = new Solver.NewmarkBeta_solver();
            slv.NewmarkBeta(M, K, h, dt, beta, N, wave, ref out_a, ref out_v, ref out_d);

            // grassshopper へのデータ出力　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            DA.SetDataList(0, Model);
            DA.SetDataList(1, out_a);
            DA.SetDataList(2, out_v);
            DA.SetDataList(3, out_d);
        }

        /// <summary>
        /// GUIDの設定
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("419c3a3a-cc48-4717-9cef-5f5647a5ecfc"); }
        }

        /// <summary>
        /// アイコンの設定。24x24 pixelsが推奨
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return _1dofResponseAnalysis.Properties.Resource.icon;
            }
        }
    }

    public class CalcNaturalPeriodComponet : GH_Component
    {
        public CalcNaturalPeriodComponet()
            : base("Calc Natural Period",                  // 名称
                   "Calc T",                               // 略称
                   "Calculation of the Natural Period",    // コンポーネントの説明
                   "rgkr",                                 // カテゴリ(タブの表示名)
                   "Response Analysis"                     // サブカテゴリ(タブ内の表示名)
                  )
        {
        }

        /// <summary>
        /// インプットパラメータの登録
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Mass", "M", "Lumped Mass(ton)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Stiffness", "K", "Spring Stiffness(kN/m)", GH_ParamAccess.item);
        }

        /// <summary>
        /// アウトプットパラメータの登録
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("NaturalPeriod", "T", "output Natural Period(sec)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Naturalfrequency", "f", "output Natural Frequency(Hz)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Naturalangular frequency", "omega", "output Natural Angular Frequency(rad/sec)", GH_ParamAccess.item);
        }
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // パラメータの定義 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            double M = double.NaN;        // 質量
            double K = double.NaN;        // 剛性
            double T = double.NaN;        // 固有周期
            double f = double.NaN;        // 固有周波数
            double omega = double.NaN;    // 固有各振動数

            // grasshopper からデータ取得　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            if (!DA.GetData(0, ref M)) { return; }
            if (!DA.GetData(1, ref K)) { return; }

            // 各値の計算 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            omega = Math.Sqrt(K / M);
            T = 2.0 * Math.PI / omega;
            f = 1.0 / T;

            // grassshopper へのデータ出力　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            DA.SetData(0, T);
            DA.SetData(1, f);
            DA.SetData(2, omega);
        }

        /// <summary>
        /// GUIDの設定
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("419c3a3a-cc48-4823-9cef-5f5647a5ecfc"); }
        }

        /// <summary>
        /// アイコンの設定。24x24 pixelsが推奨
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return _1dofResponseAnalysis.Properties.Resource.calcT_icon;
            }
        }
    }

    public class MTtoKComponet : GH_Component
    {
        public MTtoKComponet()
            : base("Calc K from M & T",                    // 名称
                   "MTtoK",                                // 略称
                   "Calculate Stiffness from Mass and Period",    // コンポーネントの説明
                   "rgkr",                                 // カテゴリ(タブの表示名)
                   "Response Analysis"                     // サブカテゴリ(タブ内の表示名)
                  )
        {
        }

        /// <summary>
        /// インプットパラメータの登録
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Mass", "M", "Lumped Mass(ton)", GH_ParamAccess.item);
            pManager.AddNumberParameter("NaturalPeriod", "T", "Natural Period(sec)", GH_ParamAccess.item);
        }

        /// <summary>
        /// アウトプットパラメータの登録
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Stiffness", "K", "Spring Stiffness(kN/m)", GH_ParamAccess.item);
        }
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // パラメータの定義 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            double M = double.NaN;        // 質量
            double K = double.NaN;        // 剛性
            double T = double.NaN;        // 固有周期

            // grasshopper からデータ取得　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            if (!DA.GetData(0, ref M)) { return; }
            if (!DA.GetData(1, ref T)) { return; }

            // 各値の計算 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            K = 4.0 * Math.PI * Math.PI / T * T * M;

            // grassshopper へのデータ出力　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            DA.SetData(0, K);
        }

        /// <summary>
        /// GUIDの設定
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("419c3a3a-cc48-4825-9cef-5f5647a5ecfc"); }
        }

        /// <summary>
        /// アイコンの設定。24x24 pixelsが推奨
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return _1dofResponseAnalysis.Properties.Resource.calcK_icon;
            }
        }
    }

    public class KTtoMComponet : GH_Component
    {
        public KTtoMComponet()
            : base("Calc M from K & T",                    // 名称
                   "KTtoM",                                // 略称
                   "Calculate Mass from Stiffness and Period",    // コンポーネントの説明
                   "rgkr",                                 // カテゴリ(タブの表示名)
                   "Response Analysis"                     // サブカテゴリ(タブ内の表示名)
                  )
        {
        }

        /// <summary>
        /// インプットパラメータの登録
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Stiffness", "K", "Spring Stiffness(kN/m)", GH_ParamAccess.item);
            pManager.AddNumberParameter("NaturalPeriod", "T", "Natural Period(sec)", GH_ParamAccess.item);
        }

        /// <summary>
        /// アウトプットパラメータの登録
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Mass", "M", "Lumped Mass(ton)", GH_ParamAccess.item);
        }
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // パラメータの定義 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            double M = double.NaN;        // 質量
            double K = double.NaN;        // 剛性
            double T = double.NaN;        // 固有周期

            // grasshopper からデータ取得　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            if (!DA.GetData(0, ref K)) { return; }
            if (!DA.GetData(1, ref T)) { return; }

            // 各値の計算 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            M = K * ( T * T ) / (4.0 * Math.PI * Math.PI);

            // grassshopper へのデータ出力　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            DA.SetData(0, M);
        }

        /// <summary>
        /// GUIDの設定
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("419c3a3a-cc48-4830-9cef-5f5647a5ecfc"); }
        }

        /// <summary>
        /// アイコンの設定。24x24 pixelsが推奨
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return _1dofResponseAnalysis.Properties.Resource.calcM_icon;
            }
        }
    }
}

/// <summary>
/// 解析関連
/// </summary>
namespace Solver
{
    /// <summary>
    /// Newmarkβ法で応答解析を行うクラス
    /// </summary>
    public class NewmarkBeta_solver
    {
        public void NewmarkBeta(double m, double k, double h, double dt, double beta, int N, double[] Ag, ref double[] out_a, ref double[] out_v, ref double[] out_d)
        {
            // 解析関連パラメータ-----------------------------
            double a = 0.0, v = 0.0, x = 0.0, an = 0.0, vn = 0.0;
            double a0 = Ag[0];                   // 初期加速度
            double v0 = 0.0;                     // 初期速度
            double d0 = 0.0;                     // 初期変位
            double c = 2 * h * Math.Sqrt(m * k); // 粘性減衰定数

            for (int n = 0; n < N; n++)
            {
                if (n == 0)  // t = 0 の時
                {
                    a = a0;
                    v = v0;
                    x = d0;
                }
                else       //  t ≠ 0 の時
                {
                    a = -(c * (v + a * dt / 2.0) + k * (x + v * dt + a * (dt * dt) * (0.5 - beta)) + m * Ag[n])
                        / (m + c * dt / 2.0 + k * (dt * dt) * beta);
                    v = v + (1.0 / 2.0) * (a + an) * dt;
                    x = x + vn * dt + beta * (a + 2.0 * an) * (dt * dt);
                }
                // 結果を出力配列に入れる。
                out_a[n] = a;
                out_v[n] = v;
                out_d[n] = x;

                an = a;
                vn = v;
            }
        }
    }
}

namespace ModelView
{
    public class ResultViewComponent : GH_Component
    {
        public ResultViewComponent()
            : base("Result View",                           // 名称
                   "Result View",                           // 略称
                   "Response Analysis of the Single dof",   // コンポーネントの説明
                   "rgkr",                                  // カテゴリ(タブの表示名)
                   "Result"                                 // サブカテゴリ(タブ内の表示名)
                  )
        {
        }

        /// <summary>
        /// インプットパラメータの登録
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Model", "Model", "Model Data", GH_ParamAccess.list);
            pManager.AddNumberParameter("Result", "Result", "Analysis Result", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Output Number", "N", "Output Result Number", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Scale", "Sc", "Scale Model", GH_ParamAccess.item, 10);
            pManager.AddNumberParameter("High", "H", "High Model", GH_ParamAccess.item, 10);
        }

        /// <summary>
        /// アウトプットパラメータの登録
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Surface", "Srf", "Output Model Surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// 解析を実行する箇所
        /// </summary>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // パラメータの定義 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            List<double> Model = new List<double>();
            List<double> Rslt = new List<double>();
            int N = 0;
            double Sc = double.NaN;
            double H = double.NaN;
            double M = double.NaN;

            // grasshopper からデータ取得　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            if (!DA.GetDataList(0, Model)) { return; }
            if (!DA.GetDataList(1, Rslt)) { return; }
            if (!DA.GetData(2, ref N)) { return; }

            // モデルのrhino上への出力＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            Point3d FMorigin = new Point3d(0, 0, H);
            Point3d FMp1 = new Point3d(1, 0, H);
            Point3d FMp2 = new Point3d(0, 1, H);
            Plane FMplane = new Plane(FMorigin, FMp1, FMp2);
            M = Model[0];

            DA.SetData(0, new Sphere(FMplane, Sc*M));
        }

        /// <summary>
        /// GUIDの設定
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("419c3a3a-cc48-4715-9cef-5f5647a5ecfc"); }
        }

        /// <summary>
        /// アイコンの設定。24x24 pixelsが推奨
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return _1dofResponseAnalysis.Properties.Resource.ResultViewicon;
            }
        }
    }

    public class ModelViewComponent : GH_Component
    {
        public ModelViewComponent()
            : base("Model View",                           // 名称
                   "Model View",                           // 略称
                   "Response Analysis of the Single dof",   // コンポーネントの説明
                   "rgkr",                                  // カテゴリ(タブの表示名)
                   "Result"                                 // サブカテゴリ(タブ内の表示名)
                  )
        {
        }

        /// <summary>
        /// インプットパラメータの登録
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Model", "Model", "Model Data", GH_ParamAccess.item,100);
            pManager.AddNumberParameter("Scale", "Sc", "Scale Model", GH_ParamAccess.item, 10);
            pManager.AddNumberParameter("High", "H", "High Model", GH_ParamAccess.item, 3500);
        }

        /// <summary>
        /// アウトプットパラメータの登録
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddSurfaceParameter("Surface", "Srf", "Output Model Surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// 解析を実行する箇所
        /// </summary>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // パラメータの定義 ＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            double Sc = double.NaN;
            double H = double.NaN;
            double M = double.NaN;
            double K = double.NaN;
            double Model = double.NaN;

            // grasshopper からデータ取得　＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            if (!DA.GetData(0, ref Model)) { return; }
            if (!DA.GetData(1, ref Sc)) { return; }
            if (!DA.GetData(2, ref H)) { return; }

            // 質点の作成＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            Point3d FMorigin = new Point3d(0, 0, H);
            Point3d FMp1 = new Point3d(1, 0, H);
            Point3d FMp2 = new Point3d(0, 1, H);
            Plane FMplane = new Plane(FMorigin, FMp1, FMp2);
            M = Model;
            Sphere FirstMass = new Sphere(FMplane, Sc * M);

            // バネの作成＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            Point3d FSPGorigin = new Point3d(0, 0, 0);
            Point3d FSPGp1 = new Point3d(1, 0, 0);
            Point3d FSPGp2 = new Point3d(0, 1, 0);
            Plane FSPGplane = new Plane(FSPGorigin, FSPGp1, FSPGp2);
            K = Model/10.0;
            Cylinder FirstSpring = new Cylinder(new Circle(FSPGplane, Sc * K), H-Sc*M);

            // モデルのrhino上への出力＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
            var Srf = new Surface[2];
            Srf[0] = FirstSpring.ToRevSurface();
            Srf[1] = FirstMass.ToRevSurface();

            DA.SetDataList(0, Srf);
        }

        /// <summary>
        /// GUIDの設定
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("419c3a3a-cc48-4715-9cef-5f5648a5ecfc"); }
        }

        /// <summary>
        /// アイコンの設定。24x24 pixelsが推奨
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return _1dofResponseAnalysis.Properties.Resource.ModelViewicon;
            }
        }
    }
}