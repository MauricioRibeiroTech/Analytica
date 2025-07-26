import numpy as np
from scipy.ndimage import generic_filter
from collections import defaultdict
from math import log2
import os
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.ticker import PercentFormatter


class AnalisadorBacias:
    def __init__(self, arquivo_matriz):
        """Inicializa o analisador carregando a matriz de bacias."""
        self.matriz = np.loadtxt(arquivo_matriz, delimiter='\t', dtype=np.int64)
        self._validar_matriz()

    def _validar_matriz(self):
        """Verifica se a matriz contém apenas -1 ou valores positivos."""
        if not np.all((self.matriz == -1) | (self.matriz > 0)):
            raise ValueError("A matriz deve conter apenas valores -1 ou positivos.")

    def calcular_estatisticas_basicas(self):
        """Calcula estatísticas básicas das bacias."""
        mascara_bacias = self.matriz != -1
        total_pontos = np.sum(mascara_bacias)
        valores_bacias = np.unique(self.matriz[self.matriz != -1])

        bacias = {}
        for bacia in valores_bacias:
            contagem = np.sum(self.matriz == bacia)
            proporcao = contagem / total_pontos if total_pontos > 0 else 0
            bacias[bacia] = {'contagem': contagem, 'proporcao': proporcao}

        return {'total_pontos': total_pontos, 'bacias': bacias}

    def identificar_fronteiras(self):
        """Identifica pontos de fronteira entre bacias."""
        fronteira = np.zeros_like(self.matriz, dtype=np.int64)
        for i in range(1, self.matriz.shape[0] - 1):
            for j in range(1, self.matriz.shape[1] - 1):
                if self.matriz[i, j] == -1:
                    continue
                vizinhos = self.matriz[i - 1:i + 2, j - 1:j + 2]
                vizinhos_bacias = vizinhos[vizinhos != -1]
                if len(np.unique(vizinhos_bacias)) > 1:
                    fronteira[i, j] = 1

        total_fronteiras = np.sum(fronteira)
        proporcao = total_fronteiras / np.sum(self.matriz != -1) if np.sum(self.matriz != -1) > 0 else 0
        return fronteira, total_fronteiras, proporcao

    def calcular_entropia_topologica(self, tamanho_janela=3):
        """Calcula a entropia topológica."""
        bacia = np.where(self.matriz == -1, 0, self.matriz).astype(np.int64)
        contagem_padroes = defaultdict(int)

        def _contar_padroes(regiao):
            padrao = tuple(regiao.flatten())
            contagem_padroes[padrao] += 1
            return 0.0

        generic_filter(
            input=bacia,
            function=_contar_padroes,
            size=tamanho_janela,
            mode='constant',
            cval=0,
            output=np.zeros_like(bacia)
        )

        total_padroes = sum(contagem_padroes.values())
        entropia = 0.0
        for contagem in contagem_padroes.values():
            probabilidade = contagem / total_padroes
            entropia -= probabilidade * log2(probabilidade) if probabilidade > 0 else 0
        return entropia

    def calcular_entropia_fronteiras(self):
        """Calcula a entropia das fronteiras."""
        fronteira, _, _ = self.identificar_fronteiras()
        hist = np.bincount(fronteira.flatten(), minlength=2)
        prob = hist / np.sum(hist)
        return -np.sum([p * log2(p) for p in prob if p > 0])

    def calcular_probabilidades_bacias(self):
        """Calcula e retorna as probabilidades de cada bacia."""
        estatisticas = self.calcular_estatisticas_basicas()
        probabilidades = {bacia: dados['proporcao'] for bacia, dados in estatisticas['bacias'].items()}
        return probabilidades


def salvar_probabilidades_txt(probabilidades, arquivo_saida):
    """Salva as probabilidades em um arquivo .txt formatado."""
    with open(arquivo_saida, 'w') as f:
        for bacia, prob in sorted(probabilidades.items()):
            f.write(f"Bacia {bacia}: {prob:.6f}\n")


def processar_lote(diretorio, arquivo_saida='resultados_bacias.csv', arquivo_fronteiras='fronteiras.csv', arquivo_probabilidades='probabilidades_bacias.txt'):
    """Processa todos os arquivos .txt no diretório e salva os resultados em CSVs."""
    resultados = []
    dados_fronteiras = []
    todas_probabilidades = []

    # Lista e ordena os arquivos numericamente
    arquivos = [arquivo for arquivo in os.listdir(diretorio) if arquivo.endswith('.txt')]
    arquivos_ordenados = sorted(arquivos, key=lambda x: int(x.split('_')[0]))

    for arquivo in arquivos_ordenados:
        caminho = os.path.join(diretorio, arquivo)
        try:
            analisador = AnalisadorBacias(caminho)
            estatisticas = analisador.calcular_estatisticas_basicas()
            entropia_top = analisador.calcular_entropia_topologica()
            entropia_front = analisador.calcular_entropia_fronteiras()
            fronteira, total_front, prop_front = analisador.identificar_fronteiras()
            
            # Calcula probabilidades para este arquivo
            probabilidades = analisador.calcular_probabilidades_bacias()
            todas_probabilidades.append((arquivo, probabilidades))

            # Dados para o CSV principal
            bacias_str = '|'.join([f"{k}:{v['contagem']}" for k, v in estatisticas['bacias'].items()])
            proporcoes_str = '|'.join([f"{k}:{v['proporcao']:.4f}" for k, v in estatisticas['bacias'].items()])

            resultados.append({
                'Arquivo': arquivo,
                'Total_Pontos': estatisticas['total_pontos'],
                'Bacias': bacias_str,
                'Proporcoes': proporcoes_str,
                'Fronteiras': total_front,
                'Proporcao_Fronteiras': prop_front,
                'Entropia_Topologica': entropia_top,
                'Entropia_Fronteiras': entropia_front
            })

            # Dados específicos para o CSV de fronteiras
            dados_fronteiras.append({
                'Arquivo': arquivo,
                'Total_Pontos': estatisticas['total_pontos'],
                'Total_Fronteiras': total_front,
                'Proporcao_Fronteiras': prop_front,
                'Densidade_Fronteiras': total_front / estatisticas['total_pontos'] if estatisticas['total_pontos'] > 0 else 0
            })

        except Exception as e:
            print(f"Erro ao processar {arquivo}: {e}")

    # Salva CSV principal
    if resultados:
        with open(arquivo_saida, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=resultados[0].keys())
            writer.writeheader()
            writer.writerows(resultados)
        print(f"Resultados gerais salvos em {arquivo_saida}")

    # Salva CSV específico de fronteiras
    if dados_fronteiras:
        with open(arquivo_fronteiras, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=dados_fronteiras[0].keys())
            writer.writeheader()
            writer.writerows(dados_fronteiras)
        print(f"Dados de fronteiras salvos em {arquivo_fronteiras}")

    # Salva probabilidades em arquivo TXT
    if todas_probabilidades:
        with open(arquivo_probabilidades, 'w') as f:
            for arquivo, probs in todas_probabilidades:
                f.write(f"=== Arquivo: {arquivo} ===\n")
                for bacia, prob in sorted(probs.items()):
                    f.write(f"Bacia {bacia}: {prob:.6f}\n")
                f.write("\n")  # Espaço entre arquivos
        print(f"Probabilidades das bacias salvas em {arquivo_probabilidades}")


def gerar_grafico_fronteiras(arquivo_csv='fronteiras.csv'):
    """Gera gráfico de proporção de fronteiras mais atrativo."""
    try:
        # Verifica estilos disponíveis e usa um alternativo se 'seaborn' não estiver disponível
        available_styles = plt.style.available
        preferred_styles = ['seaborn-v0_8', 'seaborn', 'ggplot', 'seaborn-whitegrid']

        for style in preferred_styles:
            if style in available_styles:
                plt.style.use(style)
                break
        else:
            plt.style.use('classic')  # fallback para o estilo clássico

        df = pd.read_csv(arquivo_csv)
        fig, ax = plt.subplots(figsize=(14, 7))

        # Gráfico de barras com cores gradientes
        colors = plt.cm.viridis(df['Proporcao_Fronteiras'] * 1.2)
        bars = ax.bar(df['Arquivo'], df['Proporcao_Fronteiras'],
                      color=colors,
                      edgecolor='black', linewidth=0.7)

        # Adiciona rótulos de valor
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height,
                    f'{height:.2%}', ha='center', va='bottom',
                    fontsize=9)

        # Formatação do eixo Y em porcentagem
        ax.yaxis.set_major_formatter(PercentFormatter(1.0))

        # Títulos e labels
        ax.set_title('Proporção de Pontos de Fronteira por Arquivo',
                     fontsize=16, pad=20)
        ax.set_xlabel('Arquivo', fontsize=12)
        ax.set_ylabel('Proporção de Fronteiras', fontsize=12)

        # Rotação dos rótulos do eixo X
        plt.xticks(rotation=45, ha='right')

        # Grid e ajustes finais
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()

        # Salva em alta resolução
        plt.savefig('proporcao_fronteiras_visual.png', dpi=300, bbox_inches='tight')
        plt.close()

        print("Gráfico de fronteiras salvo como 'proporcao_fronteiras_visual.png'")

    except Exception as e:
        print(f"Erro ao gerar gráfico de fronteiras: {e}")
        raise


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Uso: python analisador_bacias_lote.py <diretorio> [saida_geral.csv] [saida_fronteiras.csv] [saida_probabilidades.txt]")
        sys.exit(1)

    diretorio = sys.argv[1]
    arquivo_saida = sys.argv[2] if len(sys.argv) > 2 else 'resultados_bacias.csv'
    arquivo_fronteiras = sys.argv[3] if len(sys.argv) > 3 else 'fronteiras.csv'
    arquivo_probabilidades = sys.argv[4] if len(sys.argv) > 4 else 'probabilidades_bacias.txt'

    processar_lote(diretorio, arquivo_saida, arquivo_fronteiras, arquivo_probabilidades)

    if os.path.exists(arquivo_fronteiras):
        gerar_grafico_fronteiras(arquivo_fronteiras)