# Free Open-Source Windows Code Signing (SignPath Foundation)

This guide documents the automated code-signing setup using the **[SignPath Foundation](https://signpath.org/)**, which provides **100% free** cloud HSM Authenticode code signing for qualifying open-source software projects.

---

## 1. Why SignPath?

* **Cost:** **$0.00 (Zero)** — Funded by the SignPath Foundation for open-source software.
* **No Subscriptions / Hardware Tokens:** Keys are protected in cloud HSMs; no Azure subscriptions, physical USB tokens, or third-party certificate fees.
* **Eliminates SmartScreen:** Generates legitimate Authenticode signatures backed by recognized Certificate Authorities (e.g., Sectigo / DigiCert root CAs), removing the blue *"Windows protected your PC / Unknown Publisher"* dialog.

---

## 2. One-Time Onboarding Steps

### Step 1: Submit Project Application
1. Visit the [SignPath Foundation](https://signpath.org/) and click **Apply for free code signing**.
2. Provide your repository details:
   * **Project Name:** Correlation
   * **Repository URL:** `https://github.com/Isurwars/Correlation`
   * **License:** GNU Affero General Public License v3.0 (`AGPL-3.0`)
3. Once approved (typically 1-3 business days), you will receive access to your SignPath organization dashboard.

### Step 2: Configure Project and Policy in SignPath
1. In the SignPath dashboard, navigate to **Projects** > **Add Project**.
   * **Project Slug:** `correlation`
2. Under **Artifact Configurations**, configure the target format:
   * **Artifact Type:** `pe-file` (PE executable / NSIS installer `*.exe`)
3. Create a **Signing Policy**:
   * **Signing Policy Slug:** `release-signing`
   * **Certificate Type:** Open Source Code Signing Certificate
   * **Approval:** Set to automated CI signing for tag pushes / releases.

### Step 3: Generate API Token
1. Go to **Organization Settings** > **API Tokens**.
2. Generate an API token with permissions to submit signing requests for the `correlation` project.

---

## 3. GitHub Repository Configuration

In your GitHub repository, navigate to **Settings** > **Secrets and variables** > **Actions**:

### Repository Secrets
| Secret Name | Description | Example |
| :--- | :--- | :--- |
| `SIGNPATH_API_TOKEN` | API token generated in SignPath dashboard | `sp_live_...` |
| `SIGNPATH_ORGANIZATION_ID` | Your SignPath Organization GUID / ID | `12345678-abcd-ef01-2345-6789abcdef01` |

### Optional Repository Variables
| Variable Name | Description | Default |
| :--- | :--- | :--- |
| `SIGNPATH_PROJECT_SLUG` | Project slug configured on SignPath | `correlation` |
| `SIGNPATH_SIGNING_POLICY_SLUG` | Signing policy slug configured on SignPath | `release-signing` |

---

## 4. Release Workflow Automation

The release pipeline ([`.github/workflows/release.yml`](../.github/workflows/release.yml)) handles the signing process automatically:

1. **Detection:** Step `check-signpath` checks for `SIGNPATH_API_TOKEN` and `SIGNPATH_ORGANIZATION_ID`.
2. **If Configured:**
   * Uploads `build/correlation-*-win64.exe` as a temporary GitHub workflow artifact (`unsigned-installer`).
   * Submits a signing request to SignPath via `signpath/github-action-submit-signing-request@v3`.
   * Waits for SignPath cloud HSM signing and downloads the signed `.exe` to `build/signed/`.
   * Replaces the release binary with the signed version.
3. **If Pending / Not Configured:**
   * Emits a GitHub Actions `::warning::` notice and releases the unsigned installer without failing CI.

